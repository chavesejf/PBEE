#!/bin/python3

# ===================================================================================
# Name.......: PBEE - Protein Binding Energy Estimator
# Authors....: Roberto D. Lins, Elton J. F. Chaves, and João Sartori
# Contact....: linsrd@gmail.com
# Description: A pipeline that used ML model based on Rosetta descriptors to predict
#              the binding affinity of protein-protein complexes.
# ===================================================================================
from modules.detect_ions  import *
from modules.detect_gaps  import *
from modules.superlearner import *
from modules.rosetta_descriptors import Get_descriptors
from shutil import which
from pyrosetta import *
import pandas as pd
import os
import math
import time
import glob
import shutil
import argparse
import subprocess

def pre_processing(pdbfiles):
    bad_structures = []
    for mol, pdb in enumerate(pdbfiles):
        basename = os.path.basename(pdb[:-4])
        outdir = f'{args.odir[0]}/pbee_outputs/{basename}'

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
                print_infos(message=f'[{mol}] warning: {gap} gap(s) found.', type='info')
                total_gaps += gap
        if total_gaps > 0 and frcmod_struct is False:
            bad_structures.append(pdb)
            shutil.rmtree(outdir); continue
    return bad_structures

def post_processing(pdbfiles, partner1, partner2, trainedmodels, mlmodel, st):
    for mol, pdb in enumerate(pdbfiles):
        basename = os.path.basename(pdb[:-4])
        outdir = f'{args.odir[0]}/pbee_outputs/{basename}'

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
                for i, line in enumerate(lines):
                    if line.startswith('ATOM') and line[21] == ion[1][21]:
                        index = i
                lines.insert(index + 1, ion[1])
            with open(_pdb, 'w') as f:
                f.writelines(lines)

        # 3. executa o protocolo scorejd2
        _pdb = scorejd2(_pdb, basename, outdir)
        
        # 4. previne erros no rosetta
        _pdb, total_atoms = preventing_errors(_pdb, basename, outdir)

        # 6. executa o protocolo de minimização e calcula descritores de interface
        # ------------------------------------------------------------------------
        train_file_columns = pd.read_csv(f"{PbeePATH}/trainedmodels/{version}/{version}__pbee_train_file.csv")
        train_file_columns = train_file_columns.drop(columns=['pdb', 'database', 'partner1', 'partner2', 'dG_exp'])
        columns_to_remove  = ['pdb', 'database', 'partner1', 'partner2', 'dG_exp']
        x_train = pd.read_csv(f'{PbeePATH}/trainedmodels/{version}/{version}__pbee_train_file.csv', delimiter=',').drop(columns=columns_to_remove)
        y_train = pd.read_csv(f'{PbeePATH}/trainedmodels/{version}/{version}__pbee_train_file.csv', delimiter=',')['dG_exp']
        
        # ---
        if not os.path.isfile(f'{outdir}/dG_pred.csv'):
            print_infos(message=f'[{mol}] geometry optimization and interface analysis', type='protocol')
            pose, rosetta_features = Get_descriptors(_pdb, ions, outdir, basename, partner1, partner2)
            selected_columns = [col for col in train_file_columns if col in rosetta_features.columns]
            rosetta_features = rosetta_features[selected_columns]

            # salva arquivo .pdb
            if len(ions) != 0:
                pose.dump_pdb(f'{outdir}/{basename}_ions_rlx.pdb')
            else:
                pose.dump_pdb(f'{outdir}/{basename}_rlx.pdb')
            
            # checkpoint
            condition = not rosetta_features.applymap(lambda x: '-nan' in str(x)).any().any()
            if condition is False or rosetta_features['ifa_sc_value'][0] == -1:
                print_infos(message=f'[{mol}] an incorrect descriptor was found, ignoring the structure to avoid errors', type='protocol')
                continue

            # -------------
            # 7. calcula dG
            # -------------
            if frcmod_scores is False:
                outliers = detect_outliers(x_train, rosetta_features, mol)
                if outliers != 0:
                    continue
        else:
            rosetta_features = pd.read_csv(f'{outdir}/dG_pred.csv', delimiter=',')
            selected_columns = [col for col in train_file_columns if col in rosetta_features.columns]
            rosetta_features = rosetta_features[selected_columns]

        # ---
        print_infos(message=f'[{mol}] calculating ΔG[bind]', type='protocol')
        dG_pred = predictor(trainedmodels, mlengine, mlmodel, x_train, y_train, rosetta_features, columns_to_remove)
        affinity = calc_affinity(dG_pred)
        print_dG(mol, dG_pred, affinity)
        total_time = processing_time(st)

        # ---
        rosetta_features.insert(0, 'pdb',             basename)
        rosetta_features.insert(1, 'dG_pred',         dG_pred)
        rosetta_features.insert(2, 'affinity',        affinity)
        rosetta_features.insert(3, 'mlengine',        mlengine)
        rosetta_features.insert(4, 'total_atoms',     total_atoms)
        rosetta_features.insert(5, 'processing_time', total_time)
        rosetta_features.to_csv(f'{outdir}/dG_pred.csv', index=False)

        # 8 Apaga arquivos temporários
        remove_files(files=[
            glob.glob(f'{outdir}/*fasta'),
            f'{outdir}/{basename}_jd2_01.pdb',
            f'{outdir}/{basename}_jd2_02.pdb'])

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
    # ---
    atoms = 0
    pdb2  = f'{outdir}/{basename}_jd2_02.pdb'
    with open(pdb1, 'r') as inpfile, open(pdb2, 'w') as outfile:
        for line in inpfile:
            if not line.startswith('SSBOND'):
                outfile.write(line)
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    atoms += 1
    return pdb2, atoms

def pdbcleaner(pdbfile, basename, outdir, submit_dir, partner1, partner2):
    commands = [
        f'python {PbeePATH}/modules/clean_pdb.py {pdbfile} {partner1}',
        f'python {PbeePATH}/modules/clean_pdb.py {pdbfile} {partner2}',
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
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')
    
    pyrosetta.init(extra_options="\
    -corrections::beta_nov16 true \
    -mute core \
    -mute basic \
    -ignore_unrecognized_res \
    -output_pose_energies_table false \
    -renumber_pdb")
    pose = rosetta.core.import_pose.pose_from_file(pdbfile)
    pose.dump_pdb(f'{outdir}/{basename}_jd2_0001.pdb')
    
    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    return f'{outdir}/{basename}_jd2_0001.pdb'

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
    return elapsed_time

def print_infos(message, type):
    if type == 'info':
        print(f'            info: {message}')
    if type == 'structure':
        print(f'       structure: {message}')
    if type == 'protocol':
        print(f'        protocol: {message}')
    if type == 'none':
        print(f' {message}')

def print_dG(mol, dG_pred, affinity):
    print_infos(message=f'[{mol}] ΔG[bind] = {dG_pred:.3f} kcal/mol (KD = {affinity} M)', type='protocol')

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

def predictor(trainedmodels, mlengine, mlmodel, x, y, rosetta_features, columns_to_remove):
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
    print( ' =====================================================')
    print( '   Protein Engineering and Structural Genomic Group  ')    
    print( '          Oswaldo Cruz Foundation - FIOCRUZ          ')
    print( ' -----------------------------------------------------')
    print( '')
    print( ' ********* Protein Binding Energy Estimator **********')
    print( '')
    print( ' Authors: Roberto Lins, Elton Chaves, and João Sartori')
    print( '     DOI: 10.26434/chemrxiv-2023-zq1nj')
    print(f' Version: {version}')
    print( ' =====================================================')
    print( '')

if (__name__ == "__main__"):
    warnings.filterwarnings("ignore")
    
    # Versão do script
    version = 'v1.1'

    # Define o tempo de início do script
    st = time.time()

    # Define variável que armazena o diretório de submissão
    submit_dir = os.getcwd()

    # Imprime o cabeçalho na tela
    header(version)

    # Define PbeePATH
    PbeePATH = configure_PbeePATH()

    # Define modelo ML
    trainedmodels = configure_mlmodels(PbeePATH)
    mlmodels = {
        'sl': f'{PbeePATH}/trainedmodels/{version}/{version}__SuperLearner.pkl',
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
    parser.add_argument('--ion_dist_cutoff', nargs=1, type=float, default=[2], metavar='',
    help='int | cutoff distance (Å) to detect ion(s) close to the protein atoms (default=2)')
    parser.add_argument('--frcmod_struct', action='store_true',
    help='ignores warning messages about structure(s) with gap(s)')
    parser.add_argument('--frcmod_scores', action='store_true',
    help='ignores warning messages about low-quality descriptors')

    # ---
    args            = parser.parse_args()
    pdbfiles        = args.ipdb
    partner1        = args.partner1[0]
    partner2        = args.partner2[0]
    odir            = args.odir[0]
    mlengine        = args.mlengine[0]
    mlmodel         = mlmodels[mlengine]
    ion_dist_cutoff = args.ion_dist_cutoff[0]
    frcmod_struct   = args.frcmod_struct
    frcmod_scores   = args.frcmod_scores
    
    # Mostra parâmetros do script na tela
    # -----------------------------------
    print(f'        mlengine: {mlmodel}')
    print(f'      output_dir: {odir}')
    print(f'        partner1: {partner1}')
    print(f'        partner2: {partner2}')
    print(f' ion_dist_cutoff: {ion_dist_cutoff}')
    if frcmod_struct is True:
        print(f'   frcmod_struct: {frcmod_struct}')
    if frcmod_scores is True:
        print(f'   frcmod_scores: {frcmod_scores}')

    # Pré-processamento
    # -----------------
    bad_structures = pre_processing(pdbfiles)
    if frcmod_struct is False:
        pdbfiles = [item for item in pdbfiles if item not in bad_structures]
    else:
        pass
    
    # Pós-processamento
    # -----------------
    print_infos(message=f'total structures: {len(pdbfiles)}', type='info')
    if len(pdbfiles) != 0:
        post_processing(pdbfiles, partner1, partner2, trainedmodels, mlmodel, st)
    else:
        print_infos(message='nothing to do', type='info'); print_end()

    elapsed_time = processing_time(st)
    print(' processing time:', time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
    print_end()
