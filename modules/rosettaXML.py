#!/bin/python3
# Written by Elton Chaves

import os, string

def prepareXML(outdir, prefix, partner1, partner2):
    partner1_index = []
    if len(partner1) != 1:
        partner1_chains = []
        for index, item in enumerate(partner1):
            jump = index + 1
            partner1_index.append(jump)
            chain = string.ascii_uppercase[index]
            partner1_chains.append(chain)
        partner1_index_str = ','.join(map(str, partner1_index))
        partner1 = ','.join(map(str, partner1_chains))
    else:
        jump = 1
        partner1_index_str = '1'
        partner1 = string.ascii_uppercase[0]
    
    partner2_index = []
    if len(partner2) != 1:
        partner2_chains = []
        for index, item in enumerate(partner2):
            partner2_index.append(index + jump + 1)
            chain = string.ascii_uppercase[index + jump]
            partner2_chains.append(chain)
        partner2_index_str = ','.join(map(str, partner2_index))
        partner2 = ','.join(map(str, partner2_chains))
    else:
        partner2_index_str = jump + 1
        partner2 = string.ascii_uppercase[jump]
    interface = partner1.replace(',','') + '_' + partner2.replace(',','')

    for entry in os.listdir(outdir):
        if entry.endswith('.pdb'):
            pdb = entry
    xml = f'{outdir}/{prefix}.xml'
    
    # score function
    scrfxn = '<ScoreFunction name="beta" weights="beta_nov16"/>'
    
    # task operations
    task_operations = {
        'init': '<InitializeFromCommandline name="init"/>',
        'rtiv': f'<RestrictToInterfaceVector name="rtiv" chain1_num="{partner1_index_str}" chain2_num="{partner2_index_str}"/>'}
    
    # simple metrics
    simple_metrics = {
        'int_energy': '<InteractionEnergyMetric name="ie" residue_selector="partner1" residue_selector2="partner2" scorefxn="beta"/>'}

    # filters
    filters = {
        'sc': '<ShapeComplementarity name="sc" residue_selector1="partner1" residue_selector2="partner2" confidence="0"/>',
        'cms': '<ContactMolecularSurface name="cms" distance_weight="0.5" target_selector="partner1" binder_selector="partner2" confidence="0"/>',
        'holes': f'<InterfaceHoles name="holes" jump="{jump}" confidence="0"/>',
        'ddg_filter1': f'<Ddg name="ddg_filter1" scorefxn="beta" jump="{jump}" chain_num="{partner2_index_str}" repeats="1" repack="0" repack_bound="0" repack_unbound="0" threshold="99999" confidence="0"/>'}
    
    # movers
    movers = {
        'min_sdchains': f'<MinMover name="min1" scorefxn="beta" jump="{jump}" max_iter="50000" tolerance="0.0001" cartesian="0" bb="0" chi="1" bb_task_operations="init" chi_task_operations="init"/>',
        'min_sdchains_int': f'<MinMover name="min1" scorefxn="beta" jump="{jump}" max_iter="50000" tolerance="0.0001" cartesian="0" bb="0" chi="1" bb_task_operations="init" chi_task_operations="init,rtiv"/>',
        'min_bb': f'<MinMover name="min2" scorefxn="beta" jump="{jump}" max_iter="50000" tolerance="0.0001" cartesian="0" bb="1" chi="1" bb_task_operations="init" chi_task_operations="init"/>',
        'min_bb_int': f'<MinMover name="min2" scorefxn="beta" jump="{jump}" max_iter="50000" tolerance="0.0001" cartesian="0" bb="1" chi="1" bb_task_operations="init" chi_task_operations="init,rtiv"/>',
        'ifa': f'<InterfaceAnalyzerMover name="ifa" scorefxn="beta" interface="{interface}" packstat="1" interface_sc="1" tracer="1" scorefile_reporting_prefix="ifa"/>',
        'int_energy': '<RunSimpleMetrics name="iesum" metrics="ie"/>'}
    
    # protocols
    protocol_movers = {
        'min1': '<Add mover_name="min1"/>',
        'min2': '<Add mover_name="min2"/>',
        'ifa': '<Add mover_name="ifa"/>',
        'int_energy': '<Add mover_name="iesum"/>'}
    protocol_filters = {
        'cms': '<Add filter_name="cms"/>',
        'sc': '<Add filter_name="sc"/>',
        'ddg_filter1': '<Add filter_name="ddg_filter1"/>',
        'ddg_filter2': '<Add filter_name="ddg_filter2"/>',
        'holes': '<Add filter_name="holes"/>'}

    with open(xml, 'w') as f:
        f.write('<ROSETTASCRIPTS>\n')
        f.write('   <SCOREFXNS>\n')
        f.write(f'      {scrfxn}\n')
        f.write('   </SCOREFXNS>\n')
        f.write('   <RESIDUE_SELECTORS>\n')
        f.write(f'      <Chain name="partner1" chains="{partner1}"/>\n')
        f.write(f'      <Chain name="partner2" chains="{partner2}"/>\n')
        f.write('   </RESIDUE_SELECTORS>\n')
        f.write('   <TASKOPERATIONS>\n')
        f.write(f'      {task_operations["init"]}\n')
        f.write('   </TASKOPERATIONS>\n')
        f.write('   <SIMPLE_METRICS>\n')
        f.write(f'      {simple_metrics["int_energy"]}\n')
        f.write('   </SIMPLE_METRICS>\n')
        f.write('   <FILTERS>\n')
        f.write(f'      {filters["sc"]}\n')
        f.write(f'      {filters["cms"]}\n')
        f.write(f'      {filters["holes"]}\n')
        f.write(f'      {filters["ddg_filter1"]}\n')
        f.write('   </FILTERS>\n')
        f.write('   <MOVERS>\n')
        f.write(f'      {movers["min_sdchains"]}\n')
        f.write(f'      {movers["min_bb"]}\n')
        f.write(f'      {movers["ifa"]}\n')
        f.write(f'      {movers["int_energy"]}\n')
        f.write('   </MOVERS>\n')
        f.write('   <PROTOCOLS>\n')
        f.write(f'      {protocol_movers["min1"]}\n')
        f.write(f'      {protocol_movers["min2"]}\n')
        f.write(f'      {protocol_movers["int_energy"]}\n')
        f.write(f'      {protocol_movers["ifa"]}\n')
        f.write(f'      {protocol_filters["ddg_filter1"]}\n')
        f.write(f'      {protocol_filters["cms"]}\n')
        f.write(f'      {protocol_filters["sc"]}\n')
        f.write(f'      {protocol_filters["holes"]}\n')
        f.write('   </PROTOCOLS>\n')
        f.write('</ROSETTASCRIPTS>')
    return xml
