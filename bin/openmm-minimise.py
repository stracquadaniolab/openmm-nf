#!/usr/bin/env python3
"""openmm-minimise

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmm-minimise.py [--i=<pdb>] [--max-iter=<iter>] [--tol==<tol>]
Options:
--i=<pdb>          Input the fixed PDB files from the Mutant Maker process
--max-iter=<iter>  Maximum number of iterations to arrive at minimised energy
--tol=<tol>        Tolerance level of energy change after which the protein's structure is considered minimised

"""
import logging
from docopt import docopt
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
from sys import stdout

def get_pdb(filename: str):
    pdb = PDBFile(filename)
    return pdb

def setup_forcefield():
    forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
    return forcefield

def setup_system(pdb, forcefield):
    system = forcefield.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    return system

def setup_simulation(pdb, system):
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    return simulation

def main():
    arguments = docopt(__doc__, version='openmm-minimise.py')
    pdb = get_pdb(arguments['--i'])
    forcefield = setup_forcefield()
    system = setup_system(pdb, forcefield, arguments['--no-restraints'])
    str_name = arguments['--i']
    stem = str_name.replace("_fixed.pdb","")
    #csvunfolded = str(stem + "_unfolded.csv")
    pdbunfolded = str(stem + "_unfolded.pdb")
    pdbfolded = str(stem + "_folded.pdb")
    simulation = setup_simulation(pdb, system)
    init_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
            % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    simulation.minimizeEnergy(tolerance=0.01)
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    logging.info("Minimised potential energy = %.9f kcal/mol"
            % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    PDBFile.writeFile(
        simulation.topology, final_state.getPositions(), open(pdbunfolded, "w"))
    pdbII = get_pdb(pdbunfolded)
    forcefield = setup_forcefield()
    systemII = setup_system(pdbII, forcefield)
    simulationII = setup_simulation(pdbII, systemII)
    init_state = simulationII.context.getState(getEnergy=True, getPositions=True)
    logging.info("Starting potential energy = %.9f kcal/mol"
            % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    simulationII.minimizeEnergy(tolerance=0.00000001)
    final_state = simulationII.context.getState(getEnergy=True, getPositions=True)
    logging.info("Minimised potential energy = %.9f kcal/mol"
            % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
    PDBFile.writeFile(
        simulationII.topology, final_state.getPositions(), open(pdbfolded, "w"))
if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmm-minimise.py')
    logging.basicConfig(stream=stem + ".out", level=logging.INFO)
    main()
