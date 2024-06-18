#!/bin/python3
# Written by Elton Chaves

import math

def detect_ions(pdbfile,cutoff,chains):
    ion_types = ['MG','CA','NA','CL','FE','K','ZN','MN']
    ptn_coords = []
    ion_coords_all = []
    ion_coords_str = []
    with open(pdbfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('ATOM'):
                ptn_chain = line[20:22].strip()
                if ptn_chain in chains:
                    ptn_coords.append([ptn_chain,line])
            if line.__contains__('HETATM'):
                for item in ion_types:
                    if line[70:80].__contains__(item):
                        ion_chain = line[20:22].strip()
                        if ion_chain in chains:
                            ion_coords_all.append([ion_chain,line])
        for atom in ptn_coords:
            atom_x = float(atom[1][31:38])
            atom_y = float(atom[1][39:46])
            atom_z = float(atom[1][47:54])    
            for ion in ion_coords_all:
                ion_x = float(ion[1][31:38])
                ion_y = float(ion[1][39:46])
                ion_z = float(ion[1][47:54])
                # calcula distância entre o ion e átomo da proteína
                dist = calculate_dist(atom_x,atom_y,atom_z,ion_x,ion_y,ion_z)
                if dist <= cutoff:
                    ion_coords_str.append(ion)
    return ion_coords_str

def calculate_dist(x1,y1,z1,x2,y2,z2):
    dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
    return dist
