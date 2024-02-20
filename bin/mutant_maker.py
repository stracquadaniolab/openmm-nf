#!/usr/bin/env python3
from openmm.app import PDBFile
import sys
from sys import stdout
from pdbfixer import PDBFixer

pdb = PDBFixer(sys.argv[1])
pdb.addMutations(mutationChain='A', mutationResidues='19', newResidues='ALA')
pdb.findMissingResidues()
pdb.findNonstandardResidues()
pdb.replaceNonstandardResidues()
pdb.removeHeterogens(True)
pdb.findMissingAtoms()
pdb.addMissingAtoms()
pdb.addMissingHydrogens(7.0)
PDBFile.writeFile(pdb.topology, pdb.positions, open('modified_protein.pdb', 'w'))