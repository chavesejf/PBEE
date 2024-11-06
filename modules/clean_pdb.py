#!/usr/bin/env python
''' This script cleans PDBs for Rosetta by removing extraneous information, converting residues names and renumbering.

It outputs both the cleaned PDB and a fasta file of the cleaned sequence.

Required parameters are the name of the PDB you want to clean, and the chain ids of the chains you want.

The PDB name may be specified with or without the .pdb file handle and may be provided as a gziped file.
If the PDB isn't found locally, the given 4 letter code will be fetched from the internet.

Chain id: only the specified chains will be extracted. You may specify more than one: "AB" gets you chain A and B,
and "C" gets you just chain C. Special notations are "nochain" to remove chain identiry from the output, and "ignorechain"
to get all the chains.

(Script written by Phil Bradley, Rhiju Das, Michael Tyka, TJ Brunette, and James Thompson from the Baker Lab. Edits done by Steven Combs, Sam Deluca, Jordan Willis and Rocco Moretti from the Meiler Lab.)
'''

# Function of this script: "clean" raw pdb file by following tasks so that rosetta modeling becomes easier

## starts residue number at 1
## translates certain residues to their cannonical amino acid equivalents
## removes unknown residues
## removes residues with 0 occupancy
## generates a fasta file
## and leaves the 1st model among many NMR models

from __future__ import print_function
import sys
import os
from sys import argv, stderr, stdout
from os import popen, system
from os.path import exists, basename
from optparse import OptionParser

# Local package imports

from amino_acids import longer_names
from amino_acids import modres

# remote host for downloading pdbs
remote_host = ''

shit_stat_insres = False
shit_stat_altpos = False
shit_stat_modres = False
shit_stat_misdns = False  # missing density!

fastaseq = {}
pdbfile = ""


def download_pdb(pdb_id, dest_dir):
    # print("downloading %s" % ( pdb_id ))
    url = 'http://www.rcsb.org/pdb/files/%s.pdb.gz' % (pdb_id.upper())
    dest = '%s/%s.pdb.gz' % (os.path.abspath(dest_dir), pdb_id)
    wget_cmd = 'wget --quiet %s -O %s' % (url, dest)
    print( wget_cmd )
    if remote_host:
        wget_cmd = 'ssh %s %s' % (remote_host, wget_cmd)

    lines = popen(wget_cmd).readlines()
    if (exists(dest)):
        return(dest)
    else:
        print( "Error: didn't download file!" )


def check_and_print_pdb(count, residue_buffer, residue_letter):
    global pdbfile
  # Check that CA, N and C are present!def check_and_print_pdb( outid, residue_buffer )
    hasCA = False
    hasN = False
    hasC = False
    for line in residue_buffer:
        atomname = line[12:16]
        # Only add bb atoms if they have occupancy!
        occupancy = float(line[55:60])
        if atomname == " CA " and occupancy > 0.0:
            hasCA = True
        if atomname == " N  " and occupancy > 0.0:
            hasN = True
        if atomname == " C  " and occupancy > 0.0:
            hasC = True

  # if all three backbone atoms are present withoccupancy proceed to print the residue
    if hasCA and hasN and hasC:
        for line in residue_buffer:
            # add linear residue count
            newnum = '%4d ' % count
            line_edit = line[0:22] + newnum + line[27:]
            # write the residue line
            pdbfile = pdbfile + line_edit

    # finally print residue letter into fasta strea
        chain = line[21]
        try:
            fastaseq[chain] += residue_letter
        except KeyError:
            fastaseq[chain] = residue_letter
    # count up residue number
        count = count + 1
        return True
    return False

def get_pdb_filename( name ):
    '''Tries various things to get the filename to use.
    Returns None if no acceptable file exists.'''
    if( os.path.exists( name ) ):
        return name
    if( os.path.exists( name + '.pdb' ) ):
        return name + '.pdb'
    if( os.path.exists( name + '.pdb.gz' ) ):
        return name + '.pdb.gz'
    if( os.path.exists( name + '.pdb1.gz' ) ):
        return name + '.pdb1.gz'
    name = name.upper()
    if( os.path.exists( name ) ):
        return name
    if( os.path.exists( name + '.pdb' ) ):
        return name + '.pdb'
    if( os.path.exists( name + '.pdb.gz' ) ):
        return name + '.pdb.gz'
    if( os.path.exists( name + '.pdb1.gz' ) ):
        return name + '.pdb1.gz'
    # No acceptable file found
    return None


def open_pdb( name ):
    '''Open the PDB given in the filename (or equivalent).
    If the file is not found, then try downloading it from the internet.

    Returns: (lines, filename_stem)
    '''
    filename = get_pdb_filename( name )
    if filename is not None:
        print( "Found existing PDB file at", filename )
    else:
        print( "File for %s doesn't exist, downloading from internet." % (name) )
        filename = download_pdb(name[0:4].upper(), '.')
        global files_to_unlink
        files_to_unlink.append(filename)

    stem = os.path.basename(filename)
    if stem[-3:] == '.gz':
        stem = stem[:-3]
    if stem[-5:] == '.pdb1':
        stem = stem[:-5]
    if stem[-4:] == '.pdb':
        stem = stem[:-4]

    if filename[-3:] == '.gz':
        lines = popen('zcat '+filename, 'r').readlines()
    else:
        lines = open(filename, 'r').readlines()

    return lines, stem

#############################################
# Program Start
#############################################

parser = OptionParser(usage="%prog [options] <pdb> <chain id>",
        description=__doc__)
parser.add_option("--nopdbout", action="store_true",
        help="Don't output a PDB.")
parser.add_option("--allchains", action="store_true",
        help="Use all the chains from the input PDB.")
parser.add_option("--removechain", action="store_true",
        help="Remove chain information from output PDB.")
parser.add_option("--keepzeroocc", action="store_true",
        help="Keep zero occupancy atoms in output.")

options, args = parser.parse_args()

if 'nopdbout' in args:
    options.nopdbout = True
    args.remove('nopdbout')

if 'ignorechain' in args:
    options.allchains = True
    #Don't remove, because we're also using it as a chain designator

if 'nochain' in args:
    options.removechain = True
    options.allchains = True
    #Don't remove, because we're also using it as a chain designator

if len(args) != 2:
    parser.error("Must specify both the pdb and the chain id")

files_to_unlink = []

if args[1].strip() != "ignorechain" and args[1].strip() != "nochain":
    chainid = args[1].upper()
else:
    chainid = args[1]

lines, filename_stem = open_pdb( args[0] )

oldresnum = '   '
count = 1

residue_buffer = []
residue_letter = ''

if chainid == '_':
    chainid = ' '

for line in lines:

    if line.startswith('ENDMDL'): break  # Only take the first NMR model
    if len(line) > 21 and ( line[21] in chainid or options.allchains):
        if line[0:4] != "ATOM" and line[0:6] != 'HETATM':
            continue

        line_edit = line
        resn = line[17:20]

        # Is it a modified residue ?
        # (Looking for modified residues in both ATOM and HETATM records is deliberate)
        if resn in modres:
            # if so replace it with its canonical equivalent !
            orig_resn = resn
            resn = modres[resn]
            line_edit = 'ATOM  '+line[6:17]+ resn + line[20:]

            if orig_resn == "MSE":
                # don't count MSE as modified residues for flagging purposes (because they're so common)
                # Also, fix up the selenium atom naming
                if (line_edit[12:14] == 'SE'):
                    line_edit = line_edit[0:12]+' S'+line_edit[14:]
                if len(line_edit) > 75:
                    if (line_edit[76:78] == 'SE'):
                        line_edit = line_edit[0:76]+' S'+line_edit[78:]
            else:
                shit_stat_modres = True

        # Only process residues we know are valid.
        if resn not in longer_names:
            continue

        resnum = line_edit[22:27]

        # Is this a new residue
        if not resnum == oldresnum:
            if residue_buffer != []:  # is there a residue in the buffer ?
                if not check_and_print_pdb(count, residue_buffer, residue_letter):
                    # if unsuccessful
                    shit_stat_misdns = True
                else:
                    count = count + 1

            residue_buffer = []
            residue_letter = longer_names[resn]

        oldresnum = resnum

        insres = line[26]
        if insres != ' ':
            shit_stat_insres = True

        altpos = line[16]
        if altpos != ' ':
            shit_stat_altpos = True
            if altpos == 'A':
                line_edit = line_edit[:16]+' '+line_edit[17:]
            else:
                # Don't take the second and following alternate locations
                continue

        if options.removechain:
            line_edit = line_edit[:21]+' '+line_edit[22:]

        if options.keepzeroocc:
            line_edit = line_edit[:55] +" 1.00"+ line_edit[60:]

        residue_buffer.append(line_edit)


if residue_buffer != []: # is there a residue in the buffer ?
    if not check_and_print_pdb(count, residue_buffer, residue_letter):
        # if unsuccessful
        shit_stat_misdns = True
    else:
        count = count + 1

flag_altpos = "---"
if shit_stat_altpos:
    flag_altpos = "ALT"
flag_insres = "---"
if shit_stat_insres:
    flag_insres = "INS"
flag_modres = "---"
if shit_stat_modres:
    flag_modres = "MOD"
flag_misdns = "---"
if shit_stat_misdns:
    flag_misdns = "DNS"

nres = len("".join(list(fastaseq.values())))

flag_successful = "OK"
if nres <= 0:
    flag_successful = "BAD"

if chainid == ' ':
    chainid = '_'

print( filename_stem, "".join(chainid), "%5d" % nres, flag_altpos,  flag_insres,  flag_modres,  flag_misdns, flag_successful )

if nres > 0:
    if not options.nopdbout:
        # outfile = string.lower(pdbname[0:4]) + chainid + pdbname[4:]
        outfile = filename_stem + "_" + chainid + ".pdb"

        outid = open(outfile, 'w')
        outid.write(pdbfile)
        outid.write("TER\n")
        outid.close()

    fastaid = stdout
    if not options.allchains:
        for chain in fastaseq:
            fastaid.write('>'+filename_stem+"_"+chain+'\n')
            fastaid.write(fastaseq[chain])
            fastaid.write('\n')
            handle = open(filename_stem+"_"+"".join(chain) + ".fasta", 'w')
            handle.write('>'+filename_stem+"_"+"".join(chain)+'\n')
            handle.write(fastaseq[chain])
            handle.write('\n')
            handle.close()
    else:
        fastaseq = ["".join(list(fastaseq.values()))]
        fastaid.write('>'+filename_stem+"_"+chainid+'\n')
        fastaid.writelines(fastaseq)
        fastaid.write('\n')
        handle = open(filename_stem+"_"+chainid + ".fasta", 'w')
        handle.write('>'+filename_stem+"_"+chainid+'\n')
        handle.writelines(fastaseq)
        handle.write('\n')
        handle.close()

if len(files_to_unlink) > 0:
    for file in files_to_unlink:
        os.unlink(file)
