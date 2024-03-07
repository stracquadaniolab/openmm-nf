#!/usr/bin/env python3
"""mutant_maker

Fix PDB files, adding missing residues, atoms and hydrogens, then create mutations in desired locations   
Usage:
mutant_maker.py [--incsv=<incsv>] [--from-col=<col>] [--in-pdb=<pdb>] [--chain=<chain>] [--pH=<pH>]

Options:
--incsv=<incsv>    Read data from csv file
--from-col=<col>   Select column title containing mutations
--in-pdb=<pdb>     PDB file of protein wildtype
--chain=<chain>    Select chain upon which to perform mutation
--pH=<ph>          Set pH of the protein  
"""

from docopt import docopt
import pandas as pd
import re
from openmm.app import PDBFile
import sys
from sys import stdout
from pdbfixer import PDBFixer

def get_data(csvname: str, col: str):
    data = pd.read_csv(csvname)
    clipped_data=pd.DataFrame(data=data, columns=[col])
    divided_data = clipped_data.Mutation.str.extract(r'([a-zA-Z]+)([^a-zA-Z]+)([a-zA-Z]+)')
    divided_cleaned = divided_data.to_string(header=None, index=False)
    amino_index = {'G': 'GLY' , 'L': 'LEU', 'R': 'ARG', 'A': 'ALA', 'Y': 'TYR', 'S': 'SER', 'M': 'MET', 'I': 'ILE', 'T': 'THR', 'C': 'CYS', 'V': 'VAL', 'P': 'PRO', 'F': 'PHE', 'W': 'TRP', 'H': 'HIS', 'K': 'LYS', 'D': 'ASP', 'E': 'GLU', 'N': 'ASN', 'Q': 'GLN', ' ': '-'}
    new_rep_codes = re.sub(r"[GLRAYSMITCVPFWHKDENQ ]", lambda x: amino_index[x.group(0)], divided_cleaned)
    new_rep_cleaned = re.sub("--","-", new_rep_codes)
    return new_rep_cleaned


def clean_wildtype(pdbname: str, pH: str):
    pH_fl = float(pH)
    pdb = PDBFixer(pdbname)
    pdb.findMissingResidues()
    pdb.findNonstandardResidues()
    pdb.replaceNonstandardResidues()
    pdb.removeHeterogens(False)
    pdb.findMissingAtoms()
    pdb.addMissingAtoms()
    pdb.addMissingHydrogens(pH_fl)
    PDBFile.writeFile(pdb.topology, pdb.positions, open("wildtype_fixed.pdb", 'w'))
    

def create_mutants(pdbname: str, new_rep_cleaned, chain: str, pH: str):
    pH_fl = float(pH)
    for mutant in new_rep_cleaned.splitlines():
        mutpdb = PDBFixer(pdbname)
        mutpdb.applyMutations([mutant], chain)
        mutpdb.findMissingResidues()
        mutpdb.findNonstandardResidues()
        mutpdb.replaceNonstandardResidues()
        mutpdb.removeHeterogens(False)
        mutpdb.findMissingAtoms()
        mutpdb.addMissingAtoms()
        mutpdb.addMissingHydrogens(pH_fl)
        PDBFile.writeFile(pdb.topology, pdb.positions, open(mutant + "_fixed.pdb", 'w'))


def main():
    arguments = docopt(__doc__, version='mutant_maker.py')
    new_rep_cleaned = get_data(arguments['--incsv'], arguments['--from-col'])
    clean_wildtype(arguments['--inpdb'], arguments['--pH'])
    create_mutants(arguments['--inpdb'], new_rep_cleaned, arguments['--chain'], arguments['--pH'])

if __name__ == '__main__':
    main()
