#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  8 11:05:09 2024

@author: joao
"""

from pyrosetta import rosetta
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.kinematics import MoveMap
from rosetta.core.kinematics import FoldTree
from rosetta.core.pack.task import TaskFactory
from rosetta.core.pack.task import operation
from rosetta.core.simple_metrics import metrics
from rosetta.core.select import residue_selector as selections
from rosetta.core import select
from rosetta.core.select.movemap import *
from rosetta.protocols import minimization_packing as pack_min
from rosetta.protocols import relax as rel
from rosetta.protocols.antibody.residue_selector import CDRResidueSelector
from rosetta.protocols.antibody import *
from rosetta.protocols.loops import *
from rosetta.protocols.relax import FastRelax
from rosetta.protocols.docking import setup_foldtree
from rosetta.protocols import *
from rosetta.core.scoring.methods import EnergyMethodOptions
import pyrosetta
import pandas as pd
import os
import sys

def jd2_format(pdbfile, basename, outdir):
    pyrosetta.init(extra_options="-corrections::beta_nov16 true -ignore_unrecognized_res -output_pose_energies_table false -renumber_pdb")
    pose = rosetta.core.import_pose.pose_from_file(pdbfile)
    pose.dump_pdb(f'{outdir}/{basename}_jd2_0001.pdb')

# Define a função de minimização
def minimize(pose, scorefxn, minimizer_type):
    # Cria um MoveMap dependendo do tipo de minimização
    movemap = pyrosetta.rosetta.core.kinematics.MoveMap()
    if minimizer_type == 'minmover1':
        movemap.set_bb(False)
        movemap.set_chi(True)
        tolerance = 0.0001
    elif minimizer_type == 'minmover2':
        movemap.set_bb(True)
        movemap.set_chi(True)
        tolerance = 0.0001

    # Cria o MinMover
    min_mover = rosetta.protocols.minimization_packing.MinMover()
    min_mover.score_function(scorefxn)
    min_mover.max_iter(50000)
    min_mover.tolerance(tolerance)
    min_mover.cartesian(False)
    min_mover.movemap(movemap)

    # Aplica a minimização à pose
    min_mover.apply(pose)

def Get_energy_per_term(pose, scorefxn):
    # Get and display the individual weighted energy terms
    score_types = rosetta.core.scoring.ScoreType
    weights = scorefxn.weights()
    energy_map = pose.energies().total_energies()  # Get total energies for the pose
    energy_terms = {}
    
    # Loop through all score types
    for score_type in rosetta.core.scoring.ScoreType.__members__.values():
        weight = weights[score_type]
        if weight != 0:  # Only include terms with non-zero weights
            term_value = energy_map[score_type] * weight
            energy_terms[score_type.name] = term_value
    return energy_terms

def Interaction_energy_metric(pose, scorefxn, partner1, partner2):
    # Create the OrResidueSelector for partner1
    partner1_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner1
    for chain in partner1:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner1_selector.add_residue_selector(chain_selector)
    
    # Create the OrResidueSelector for partner2
    partner2_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner2
    for chain in partner2:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner2_selector.add_residue_selector(chain_selector)
    
    # Initialize InteractionEnergyMetric
    interaction_energy_metric = rosetta.core.simple_metrics.metrics.InteractionEnergyMetric()
    interaction_energy_metric.set_scorefunction(scorefxn)
    
    # Set both residue selectors
    interaction_energy_metric.set_residue_selectors(partner1_selector, partner2_selector)
    
    # Calculate the interaction energy metric
    interaction_energy = interaction_energy_metric.calculate(pose)
    
    # Return the result
    return interaction_energy

def Contact_molecular_surface(pose, partner1, partner2):
    # Create the OrResidueSelector for partner1
    partner1_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner1
    for chain in partner1:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner1_selector.add_residue_selector(chain_selector)
    
    # Create the OrResidueSelector for partner2
    partner2_selector = rosetta.core.select.residue_selector.OrResidueSelector()
    
    # Add ChainSelectors for each chain in partner2
    for chain in partner2:
        chain_selector = rosetta.core.select.residue_selector.ChainSelector(chain)
        partner2_selector.add_residue_selector(chain_selector)
    
    cms_filter = rosetta.protocols.simple_filters.ContactMolecularSurfaceFilter()
    cms_filter.selector1(partner1_selector)
    cms_filter.selector2(partner2_selector)
    cms_filter.distance_weight(0.5)
    cms_filter.set_user_defined_name("cms")
    cms_filter.apply(pose)
    return cms_filter.score(pose)

def Interface_analyzer_mover(pose, partner1, partner2):
    partners = f"{partner1}_{partner2}"
    ifa_mover = rosetta.protocols.analysis.InterfaceAnalyzerMover(partners)
    ifa_mover.set_use_tracer(True)
    ifa_mover.set_compute_packstat(True)
    ifa_mover.set_scorefile_reporting_prefix("ifa")
    ifa_mover.apply(pose)
    interface_data = ifa_mover.get_all_data()
    ifa_mover.add_score_info_to_pose(pose)
    score_map = pose.scores
    # Filter the elements where the key starts with "ifa"
    ifa_data = {score_name: value for score_name, value in score_map.items() if score_name.startswith("ifa")}
    return ifa_data

def Get_interface_selector(pose, partner1, partner2):
    partners  = f"{partner1}_{partner2}"
    ifa_mover = rosetta.protocols.analysis.InterfaceAnalyzerMover(partners)
    ifa_mover.set_use_tracer(True)
    ifa_mover.set_compute_packstat(True)
    ifa_mover.set_scorefile_reporting_prefix("ifa")
    ifa_mover.apply(pose)
    interface_data = ifa_mover.get_all_data()
    residues_in_interface = []
    for i in range(1, len(interface_data.interface_residues[1]) + 1):
        if interface_data.interface_residues[1][i]:
            residues_in_interface.append(i)
    residue_selector = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    residue_selector.set_index(','.join(map(str, residues_in_interface)))
    return residue_selector

def Get_descriptors(pdb, ions, outdir, basename, partner1, partner2):
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')

    if len(ions) != 0:  
        pyrosetta.init(extra_options="\
        -mute core \
        -mute basic \
        -ex1 -ex2 -ex2aro \
        -use_input_sc \
        -flip_HNQ \
        -no_optH false \
        -corrections::beta_nov16 true \
        -output_pose_energies_table false \
        -auto_setup_metals")
    else:
        pyrosetta.init(extra_options="\
        -mute core \
        -mute basic \
        -ex1 -ex2 -ex2aro \
        -use_input_sc \
        -flip_HNQ \
        -no_optH false \
        -corrections::beta_nov16 true \
        -output_pose_energies_table false")
    pose = rosetta.core.import_pose.pose_from_file(pdb)
    
    # Create the beta_nov16 score function
    scorefxn = rosetta.core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16")
    scorefxn(pose)

    minimize(pose=pose, scorefxn=scorefxn, minimizer_type="minmover1")
    minimize(pose=pose, scorefxn=scorefxn, minimizer_type="minmover2")
    scorefxn(pose)

    per_term  = Get_energy_per_term(pose, scorefxn)
    testeIE   = Interaction_energy_metric(pose = pose, scorefxn = scorefxn, partner1 = partner1, partner2 = partner2)
    testeCMS  = Contact_molecular_surface(pose = pose, partner1 = partner1, partner2 = partner2)
    testeIFA  = Interface_analyzer_mover(pose  = pose, partner1 = partner1, partner2 = partner2)
    all_terms = per_term | testeIFA

    all_terms["cms"] = testeCMS
    all_terms["interaction_energy"] = testeIE
    all_terms["total_score"] = scorefxn(pose)
    all_terms_df = pd.DataFrame([all_terms])
    all_terms_df.rename(columns={
        "ifa_dG_separated/dSASAx100":"ifa_dG_separated_dSASAx100", 
        "ifa_dG_cross/dSASAx100":"ifa_dG_cross_dSASAx100"}, inplace=True)

    sys.stdout = sys.__stdout__
    sys.stderr = sys.__stderr__
    
    return pose, all_terms_df