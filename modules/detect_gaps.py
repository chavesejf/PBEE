#!/bin/python3
# Written by João Rodrigues (pdb-tools)

import sys

def detect_gaps(pdbfile):
    # Detects gaps between residues in the PDB file.
    fmt_GAPd = " {0[1]}:{0[3]}{0[2]} < {2:7.2f}A > {1[1]}:{1[3]}{1[2]}\n"
    fmt_GAPs = " {0[1]}:{0[3]}{0[2]} < Seq. Gap > {1[1]}:{1[3]}{1[2]}\n"
    centroid = ' CA '  # respect spacing. 'CA  ' != ' CA '
    distance_threshold = 4.0 * 4.0 
    prev_at = (None, None, None, None, (None, None, None))
    model = 0
    ngaps = 0
    with open(pdbfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line.startswith('MODEL'):
                model = int(line[10:14])
            elif line.startswith('ATOM'):
                atom_name = line[12:16]
                if atom_name != centroid:
                    continue
                resn = line[17:20]
                resi = int(line[22:26])
                chain = line[21]
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                at_uid = (model, chain, resi, resn, atom_name, (x, y, z))
                if prev_at[0] == at_uid[0] and prev_at[1] == at_uid[1]:
                    d = calculate_sq_atom_distance(at_uid[5], prev_at[5])
                    if d > distance_threshold:
                        #sys.stdout.write(fmt_GAPd.format(prev_at, at_uid, d))
                        ngaps += 1
                    elif prev_at[2] + 1 != at_uid[2]:
                        #sys.stdout.write(fmt_GAPs.format(prev_at, at_uid))
                        ngaps += 1
                prev_at = at_uid
    return ngaps

def calculate_sq_atom_distance(i,j):
    # Squared euclidean distance between two 3d points
    return (i[0] - j[0]) * (i[0] - j[0]) + \
        (i[1] - j[1]) * (i[1] - j[1]) + \
        (i[2] - j[2]) * (i[2] - j[2])
