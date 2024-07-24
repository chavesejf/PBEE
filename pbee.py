#!/bin/python3

# ===================================================================================
# Name.......: PBEE - Protein Binding Energy Estimator
# Authors....: Roberto D. Lins and Elton J. F. Chaves
# Contact....: linsrd@gmail.com
# Description: A pipeline that used ML model based on Rosetta descriptors to predict
#              the binding affinity of protein-protein complexes.
# ===================================================================================
from modules.detect_ions  import *
from modules.detect_gaps  import *
from modules.rosettaXML   import *
from modules.superlearner import *
from shutil import which
import os, math, time, shutil, argparse, subprocess, glob
import pandas as pd

def pre_processing(pdbfiles):
    bad_structures = []
    for mol, pdb in enumerate(pdbfiles):
        basename = os.path.basename(pdb[:-4])
        outdir = f'{args.odir[0]}/outputs_pbee/{basename}'

        # Verifica se o(s) arquivo(s) estão no formato PDB
        condition = ispdb(pdb)
        if condition is not False:
            print_infos(message=f'[{mol}] {pdb}\r', type='structure')
            # cria diretório para armazenar outputs
            if not os.path.isdir(outdir):
                os.makedirs(outdir)
        else:
            print_infos(message=f'invalid PDB file -> {os.path.basename(pdb)}.', type='structure')
            continue

        # 1. verifica se a estrutura contém partner1 e partner2
        # -----------------------------------------------------
        chains = partner_checker(pdb, partner1, partner2)
        if chains[0] <= 1:
            print_infos(message=f'[{mol}] argument error (--partner1/--partner2): chain ID not found ({chains[1]})', type='info')
            bad_structures.append(pdb)
            shutil.rmtree(outdir)
            if len(pdbfiles) == 1:
                print_end()
            else:
                continue
        
        # 2. verifica se a estrutura contém gaps
        # --------------------------------------
        partners = pdbcleaner(pdb, basename, outdir, submit_dir, partner1, partner2)
        gaps = []
        for partner in partners:
            n_gaps = detect_gaps(partner)
            gaps.append(n_gaps)
        total_gaps = 0
        for partner, gap in zip(partners, gaps):
            if gap != 0:
                print_infos(message=f'[{mol}] warning: {gap} gap(s) found -> {os.path.basename(partner)}', type='info')
                total_gaps += gap
        if total_gaps > 0 and allow_bad_struct is False:
            bad_structures.append(pdb)
            shutil.rmtree(outdir); continue
    return bad_structures

def post_processing(pdbfiles, partner1, partner2, trainedmodels, mlmodel):
    for mol, pdb in enumerate(pdbfiles):
        basename = os.path.basename(pdb[:-4])
        outdir = f'{args.odir[0]}/outputs_pbee/{basename}'

        # 1. concatena estruturas de partner1 e partner2
        # ----------------------------------------------
        print_infos(message=f'[{mol}] {pdb}', type='protocol')
        partners = [
            f'{outdir}/{basename}_{partner1}.pdb',
            f'{outdir}/{basename}_{partner2}.pdb']
        _pdb = concat_pdbs(outdir, basename, partner1=partners[0], partner2=partners[1])
        
        # 2. verifica se a estrutura original contém íon(s)
        # se existir, recupera as coordenadas xyz do(s) íon(s) e insere na estrutura concatenada
        # --------------------------------------------------------------------------------------
        ions = detect_ions(pdb, cutoff=ion_dist_cutoff, chains=[partner1, partner2])
        print_infos(message=f'[{mol}] total number of ions: {len(ions)}', type='protocol')  
        if len(ions) != 0:
            with open(_pdb, 'r') as f:
                lines = f.readlines()
            for ion in ions:
                for i,line in enumerate(lines):
                    if line.startswith('ATOM') and line[21] == ion[1][21]:
                        index = i
                lines.insert(index + 1, ion[1])
            with open(_pdb, 'w') as f:
                f.writelines(lines)

        # 3. executa o protocolo scorejd2
        _pdb = scorejd2(_pdb, basename, outdir)
        
        # 4. previne erros no rosetta
        _pdb = preventing_errors(_pdb, basename, outdir)

        # 5. prepara o script .xml
        _xml = prepareXML(outdir, basename, partner1, partner2)

        # 6. executa o protocolo de minimização e calcula descritores de interface
        # ------------------------------------------------------------------------
        if not os.path.isfile(f'{outdir}/dG_pred.csv'):
            print_infos(message=f'[{mol}] geometry optimization and interface analysis', type='protocol')
            rosetta_features = get_interface_features(_pdb, ions, _xml, outdir, submit_dir)

            # checkpoint
            condition = True
            with open(rosetta_features, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    if line.__contains__('-nan'):
                        condition = False
            if condition is False:
                print_infos(message=f'[{mol}] an incorrect descriptor was found, ignoring the structure to avoid errors', type='protocol')
                continue
            else:
                rosetta_features = pd.read_csv(json2csv(rosetta_features, outdir), delimiter=',').iloc[:,1:]
                if rosetta_features['sc'][0] == -1:
                    print_infos(message=f'[{mol}] an incorrect descriptor was found, ignoring the structure to avoid errors', type='protocol')
                    continue

            # -------------
            # 7. calcula dG
            # -------------
            columns_to_remove = ['decoy', 'database', 'partner1', 'partner2', 'complex_type', 'affinity', 'classifier', 'dG_exp']
            x_train = pd.read_csv(f'{PbeePATH}/train_file.csv', delimiter=',').drop(columns=columns_to_remove)
            y_train = pd.read_csv(f'{PbeePATH}/train_file.csv', delimiter=',')['dG_exp']
        
            if allow_bad_scores is False:
                outliers = detect_outliers(x_train, rosetta_features, mol)
                if outliers != 0:
                    continue

            print_infos(message=f'[{mol}] calculating dGbind', type='protocol')
            dG_pred = predictor(trainedmodels, mlmodel, x_train, y_train, rosetta_features, columns_to_remove)
            affinity = calc_affinity(dG_pred)
            rosetta_features.insert(0, 'pdb',      basename)
            rosetta_features.insert(1, 'dG_pred',  dG_pred)
            rosetta_features.insert(2, 'affinity', affinity)
            rosetta_features.insert(3, 'mlengine', mlengine)
            rosetta_features.to_csv(f'{outdir}/dG_pred.csv', index=False)

            # 7.3 Mostra os dados de dG_pred na tela
            print_infos(message=f'[{mol}] dGbind = {dG_pred:.3f} kcal/mol ({affinity} M)', type='protocol')

        else:
            data = pd.read_csv(f'{outdir}/dG_pred.csv', delimiter=',')
            dG_pred = data['dG_pred'][0]
            affinity = calc_affinity(dG_pred)
            print_infos(message=f'[{mol}] dGbind = {dG_pred:.3f} kcal/mol ({affinity} M)', type='protocol')

        # 8 Apaga arquivos temporários
        remove_files(files=[
            glob.glob(f'{outdir}/*fasta'),
            f'{outdir}/{basename}_jd2_01.pdb',
            f'{outdir}/{basename}_jd2_02.pdb',
            f'{outdir}/{basename}_jd2_0001.pdb',
            f'{outdir}/score.sc', 
            f'{outdir}/score_rlx.csv'])

def remove_files(files):
    for file in files:
        if type(file) is list:
            for item in file:
                os.remove(item)
        else:
            if os.path.exists(file):
                os.remove(file)

def partner_checker(pdbfile, partner1, partner2):
    chains, partners = detect_chains(pdbfile), list(partner1 + partner2)
    count = 0
    for chain in chains:
        if chain in partners:
            count += 1
    return count, chains

def detect_chains(pdbfile):
    chains = set()
    with open(pdbfile, 'r') as file:
        for line in file:
            if line.startswith('ATOM'):
                chain_id = line[21]
                chains.add(chain_id)
    return chains

def detect_outliers(x, rosetta_features, mol):
    count = 0
    for col in x.columns:
        for index, row in rosetta_features.iterrows():
            sup = True if row[col] > x[col].mean() + x[col].std() * 4 else False
            inf = True if row[col] < x[col].mean() - x[col].std() * 4 else False
            if sup is True or inf is True:
                print_infos(message=f'{[mol]} outlier -> {col} = {row[col]}', type='protocol')
                count += 1
    return count

def calc_affinity(dG):
    T = 298.15
    R = 8.314
    dG_J = dG * 4184
    affinity = float(f'{math.exp(dG_J / (R * T)):.6e}')
    return affinity

def json2csv(json_file, outdir):
    outfile = f'{outdir}/score_rlx.csv'
    df = pd.read_json(json_file, orient='records', lines=True)
    df.rename(columns={
        "ifa_dG_separated/dSASAx100":"ifa_dG_separated_dSASAx100",
        "ifa_dG_cross/dSASAx100":"ifa_dG_cross_dSASAx100"}, inplace=True)
    df.to_csv(outfile, index=False)
    return outfile

def preventing_errors(pdbfile, basename, outdir):
    pdb1 = f'{outdir}/{basename}_jd2_01.pdb'
    with open(pdbfile, "r") as input_file, open(pdb1, "w") as output_file:
        ter_found = False
        last_ter_index = -1
        for line_index, line in enumerate(input_file):
            if line.startswith("TER"):
                ter_found = True
                last_ter_index = line_index
        if ter_found:
            input_file.seek(0)
            for line_index, line in enumerate(input_file):
                if line_index == last_ter_index:
                    output_file.write("END" + line[3:])
                else:
                    output_file.write(line)
    pdb2 = f'{outdir}/{basename}_jd2_02.pdb'
    with open(pdb1, 'r') as inpfile, open(pdb2, 'w') as outfile:
        for line in inpfile:
            if not line.startswith('SSBOND'):
                outfile.write(line)
    return pdb2

def pdbcleaner(pdbfile, basename, outdir, submit_dir, partner1, partner2):
    commands = [
        f'python $ROSETTA3_TOOLS/protein_tools/scripts/clean_pdb.py {pdbfile} {partner1}',
        f'python $ROSETTA3_TOOLS/protein_tools/scripts/clean_pdb.py {pdbfile} {partner2}',
        f'mv {submit_dir}/{basename}_{partner1}.pdb {outdir}',
        f'mv {submit_dir}/{basename}_{partner2}.pdb {outdir}',
        f'mv {submit_dir}/{basename}_*.fasta {outdir}']
    for command in commands:
        subprocess.run(command, stdout=subprocess.PIPE, shell=True)
    return f'{outdir}/{basename}_{partner1}.pdb', f'{outdir}/{basename}_{partner2}.pdb'

def concat_pdbs(outdir, basename, partner1, partner2):
    outfile = f'{outdir}/{basename}_jd2.pdb'
    with open(partner1, 'r') as f1, open(partner2, 'r') as f2, open(outfile, 'w') as output_f:
        content1 = f1.read(); output_f.write(content1)
        content2 = f2.read(); output_f.write(content2)
    return outfile

def scorejd2(pdbfile,basename,outdir):
    command = f'\
    score_jd2.default.linuxgccrelease \
    -s {pdbfile} -renumber_pdb \
    -ignore_unrecognized_res \
    -out:pdb \
    -out:path:all {outdir} \
    -output_pose_energies_table false'
    subprocess.run(command, stdout=subprocess.PIPE, shell=True)
    return f'{outdir}/{basename}_jd2_0001.pdb'

def get_interface_features(pdbfile,ions,xml,outdir,submit_dir):
    if len(ions) != 0:
        json_file = f'{outdir}/score_ion_rlx.sc'
        command = f'\
        rosetta_scripts.default.linuxgccrelease \
        -s {pdbfile} \
        -parser:protocol {xml} \
        -holes:dalphaball $ROSETTA3/external/DAlpahBall/DAlphaBall.gcc \
        -ex1 -ex2 -ex2aro \
        -auto_setup_metals \
        -beta_nov16 \
        -use_input_sc \
        -flip_HNQ \
        -no_optH false \
        -out:suffix _ion_rlx \
        -out:file:scorefile_format json \
        -out:path:score {outdir} \
        -out:path:all {outdir} \
        -output_pose_energies_table false'
    else:
        json_file = f'{outdir}/score_rlx.sc'
        command = f'\
        rosetta_scripts.default.linuxgccrelease \
        -s {pdbfile} \
        -parser:protocol {xml} \
        -holes:dalphaball $ROSETTA3/external/DAlpahBall/DAlphaBall.gcc \
        -ex1 -ex2 -ex2aro \
        -beta_nov16 \
        -use_input_sc \
        -flip_HNQ \
        -no_optH false \
        -out:suffix _rlx \
        -out:file:scorefile_format json \
        -out:path:score {outdir} \
        -out:path:all {outdir} \
        -output_pose_energies_table false'
    subprocess.run(command, stdout=subprocess.PIPE, shell=True)
    return json_file

def ispdb(pdbfile):
    with open(pdbfile, 'r') as file:
        lines = file.readlines()
        count = 0
        for line in lines:
            if line.startswith("ATOM"): 
                count += 1
        if count != 0:
            return os.path.abspath(pdbfile)
        else:
            return False

def isdir(path):
    if os.path.isdir(path):
        return os.path.abspath(path)
    else:
        print_infos(message=f'error: path not found -> {path}', type='none'); print_end()

def istool(tool):
    return which(tool) is not None

def processing_time(st):
    sum_x = 0
    for i in range(1000000):
        sum_x += i
    time.sleep(1)
    elapsed_time = time.time() - st
    print(' processing time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))

def print_infos(message, type):
    if type == 'info':
        print(f'            info: {message}')
    if type == 'structure':
        print(f'       structure: {message}')
    if type == 'protocol':
        print(f'        protocol: {message}')
    if type == 'none':
        print(f' {message}')

def print_end():
    exit('\n --- End process ---\n')

def sl_predictions(X_test, models, meta_model):
    model_predictions = {}
    for name, model in models:
        yhat = model.predict(X_test)
        model_predictions[name] = yhat
    meta_X = np.column_stack(list(model_predictions.values()))
    super_learner_preds = meta_model.predict(meta_X)
    model_predictions["sl"] = super_learner_preds
    return model_predictions

def train_base_models(x, y, models):
    trained_models = []
    for model in models:
        model.fit(x, y)
        trained_models.append((model.__class__.__name__, model))
    return trained_models

def predictor(trainedmodels, mlmodel, x, y, rosetta_features, columns_to_remove):
    with open(mlmodel, 'rb') as f:
        meta_model = joblib.load(f)
    if mlengine != 'sl':
        model_predictions = {}
        yhat = meta_model.predict(rosetta_features.values)
        model_predictions[mlengine] = yhat
    else:
        models = [joblib.load(filename) for filename in trainedmodels]
        base_models = train_base_models(x, y, models)
        model_predictions = sl_predictions(rosetta_features, base_models, meta_model)
    model_predictions = model_predictions[mlengine][0]
    return model_predictions

def configure_PbeePATH():
    PbeePATH = os.path.dirname(__file__)
    condition = os.path.isdir(PbeePATH)
    if condition is False:
        print(' error: invalid PbeePATH'); print_end()
    else:
        return PbeePATH

def configure_requirements(PbeePATH):
    RosettaPATHS = ['ROSETTA3', 'ROSETTA3_BIN', 'ROSETTA3_TOOLS']
    count = 0
    for item in RosettaPATHS:
        condition = os.environ.get(item)
        if condition is None:
            print(f' error (rosetta_path): {item} is not an environment variable.\n')
            count += 1; continue
    if count == 0:
        r1 = glob.glob(f'{os.environ["ROSETTA3_BIN"]}/score_jd2.default.*')[0]
        r2 = glob.glob(f'{os.environ["ROSETTA3_BIN"]}/rosetta_scripts.default.*')[0]
        r3 = f'{os.environ["ROSETTA3_TOOLS"]}/protein_tools/scripts/clean_pdb.py'
        requirements = [r1, r2, r3]
        for item in requirements:
            condition = istool(item)
            if condition is False:
                print(f' error: requirement not found -> {item}\n'); exit()
    else:
        exit()

def configure_mlmodels(PbeePATH):
    trainedmodels = [
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_LinearRegression.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_ElasticNet.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_SVR.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_DecisionTreeRegressor.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_KNeighborsRegressor.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_AdaBoostRegressor.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_BaggingRegressor.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_RandomForestRegressor.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_ExtraTreesRegressor.pkl',
        f'{PbeePATH}/trainedmodels/{version}/{version}__basemodel_XGBRegressor.pkl']
    for item in trainedmodels:
        if os.path.isfile(item) is True:
            continue
        else:
            print(f' requirement not found: {item}'); print_end()

    return trainedmodels

def header(version):
    print( '')
    print( ' ====================================================')
    print( '   Protein Engineering and Structural Genomic Group  ')    
    print( '          Oswaldo Cruz Foundation - FIOCRUZ          ')
    print( ' ----------------------------------------------------')
    print( '')
    print( ' ********* Protein Binding Energy Estimator *********')
    print( '')
    print( ' Authors: Roberto Lins & Elton Chaves')
    print( '     DOI: -')
    print(f' Version: {version}')
    print( ' ====================================================')
    print( '')

if (__name__ == "__main__"):
    # Versão do script
    version = 'v1.0'

    # Define o tempo de início do script
    st = time.time()

    # Define variável que armazena o diretório de submissão
    submit_dir = os.getcwd()

    # Imprime o cabeçalho na tela
    header(version)

    # Define PbeePATH
    PbeePATH = configure_PbeePATH()
    configure_requirements(PbeePATH)

    # Define modelo ML
    trainedmodels = configure_mlmodels(PbeePATH)
    mlmodels = {
        'sl': f'{PbeePATH}/trainedmodels/v{version}/v{version}__SuperLearner.pkl',
        'lr': trainedmodels[0],
        'en': trainedmodels[1],
        'sv': trainedmodels[2],
        'dt': trainedmodels[3],
        'kn': trainedmodels[4],
        'ad': trainedmodels[5],
        'bg': trainedmodels[6],
        'rf': trainedmodels[7],
        'et': trainedmodels[8],
        'xb': trainedmodels[9]
    }
    
    # Configura argumentos do script
    # ---------------------------------
    parser = argparse.ArgumentParser()
    mandatory = parser.add_argument_group('mandatory arguments')
    
    # obrigatórios
    mandatory.add_argument('--ipdb', nargs='+', type=str, required=True, metavar='',
    help='str | input file(s) in the PDB format') 
    mandatory.add_argument('--partner1', nargs=1, type=str, required=True, metavar='',
    help='str | chain ID of the binding partner (e.g.: receptor)')
    mandatory.add_argument('--partner2', nargs=1, type=str, required=True, metavar='',
    help='str | chain ID of the binding partner (e.g.: ligand)')
    
    # opcionais
    parser.add_argument('--odir', nargs=1, type=isdir, default=[submit_dir], metavar='', 
    help=f'str | output directory (default={submit_dir})')
    parser.add_argument('--mlengine', nargs=1, type=str, default=['sl'], choices=['sl','lr','en','sv','dt','kn','ad','bg','rf','et','xb'], metavar='',
    help='str | define the machine learning engine (sl, lr, en, sv, dt, kn, ad, bg, rf, et, or xb)')
    parser.add_argument('--ion_dist_cutoff', nargs=1, type=int, default=[2], metavar='',
    help='int | cutoff distance (Å) to detect ion(s) close to the protein atoms (default=2)')
    parser.add_argument('--allow_bad_struct', action='store_true',
    help='skip warning messages about gap(s) in the structure')
    parser.add_argument('--allow_bad_scores', action='store_true',
    help='skip warning messages about bad descriptors')

    # ---
    args             = parser.parse_args()
    pdbfiles         = args.ipdb
    partner1         = args.partner1[0]
    partner2         = args.partner2[0]
    odir             = args.odir[0]
    mlengine         = args.mlengine[0]
    mlmodel          = mlmodels[mlengine]
    ion_dist_cutoff  = args.ion_dist_cutoff[0]
    allow_bad_struct = args.allow_bad_struct
    allow_bad_scores = args.allow_bad_scores
    
    # Mostra parâmetros do script na tela
    # -----------------------------------
    print(f'            info: rosetta_path -> {os.environ["ROSETTA3"]}')
    print(f'            info: rosetta_path -> {os.environ["ROSETTA3_BIN"]}')
    print(f'            info: rosetta_path -> {os.environ["ROSETTA3_TOOLS"]}')
    print(f'        mlengine: {mlmodel}')
    print(f'      output_dir: {odir}')
    print(f'        partner1: {partner1}')
    print(f'        partner2: {partner2}')
    print(f' ion_dist_cutoff: {ion_dist_cutoff}')
    if allow_bad_struct is True:
        print(f'allow_bad_struct: {allow_bad_struct}')
    if allow_bad_scores is True:
        print(f'allow_bad_scores: {allow_bad_scores}')

    # Pré-processamento
    # -----------------
    bad_structures = pre_processing(pdbfiles)
    if allow_bad_struct is False:
        pdbfiles = [item for item in pdbfiles if item not in bad_structures]
    else:
        pass
    
    # Pós-processamento
    # -----------------
    print_infos(message=f'total structures: {len(pdbfiles)}', type='info')
    if len(pdbfiles) != 0:
        post_processing(pdbfiles, partner1, partner2, trainedmodels, mlmodel)
    else:
        print_infos(message='nothing to do', type='info'); print_end()

    processing_time(st); print_end()
