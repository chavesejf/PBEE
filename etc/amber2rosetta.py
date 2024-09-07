#!/bin/python3

import os
import sys
import glob

arg1 = sys.argv[1]
arg2 = sys.argv[2]

# ---
pdb_list = glob.glob(f'{arg1}/*.pdb')

for pdb in pdb_list:
    cmd1 = f"sed 's,HIE,HIS,g; s,HID,HIS,g; s,CYX,CYS,g; s,GLH,GLU,g; s,ASH,ASN,g' {pdb} > {arg2}/{os.path.basename(pdb)}"
    os.system(cmd1)
    cmd2 = f'score_jd2.default.linuxgccrelease -s {arg2}/{os.path.basename(pdb)} \
    -ignore_zero_occupancy false \
    -out:pdb \
    -out:path:all {arg2}'
    os.system(cmd2)
    
