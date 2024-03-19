#!/usr/bin/env python3#!/usr/bin/env python3
"""openmm-minimise

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmm-minimise.py [--i=<pdb>] [--no-restraints] 

Options:
--i=<pdb>          Input the fixed PDB files from the Mutant Maker process
--max-iter=<iter>  Maximum number of iterations to arrive at minimised energy
--no-restraints    Allow movement of all atoms
"""
import logging
import csv
import pandas as pd
import numpy as np
from docopt import docopt
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
from sys import stdout
import os

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
    isFile = os.path.isfile('data.csv')
    if isFile == False:
        out = {'name':[], 'min':[], 'lowest3':[], 'lowest5':[], 'lowest10':[], 'lowest20':[]}
        inidf = pd.DataFrame(out)
        inidf.to_csv('data.csv', index=False)
    df = pd.read_csv('data.csv')
    pdb = get_pdb(arguments['--i'])
    forcefield = setup_forcefield()
    system = setup_system(pdb, forcefield, arguments['--no-restraints'])
    str_name = arguments['--i']
    stem = str_name.replace("_fixed.pdb","")
    pdbunfolded = str(stem + "_unfolded.pdb")
    pdbfolded = str(stem + "_folded.pdb")
    delta_g_list = []
    for j in range(1,81):
        simulation = setup_simulation(pdb, system)
        simulation = setup_simulation(pdb, system)
        init_state = simulation.context.getState(getEnergy=True, getPositions=True)
        logging.info("Starting potential energy = %.9f kcal/mol"
            % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
        #simulation.minimizeEnergy(tolerance=0.0000000001)
        simulation.minimizeEnergy(tolerance=2.5)
        final_state = simulation.context.getState(getEnergy=True, getPositions=True)
        logging.info("Minimised potential energy = %.9f kcal/mol"
            % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
        #final_energy = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
        PDBFile.writeFile(
        simulation.topology, final_state.getPositions(), open(pdbunfolded, "w"))
        for i in range(1,2):
            pdbII = get_pdb(pdbunfolded)
            forcefield = setup_forcefield()
            systemII = setup_system(pdbII, forcefield)
            simulationII = setup_simulation(pdbII, systemII)
            init_state = simulationII.context.getState(getEnergy=True, getPositions=True)
            init_pe = init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
            logging.info("Starting potential energy = %.9f kcal/mol"
                % init_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
            simulationII.minimizeEnergy(tolerance=0.0000000000001)
            final_state = simulationII.context.getState(getEnergy=True, getPositions=True)
            final_pe = final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)
            delta_G = init_pe - final_pe
            logging.info("Minimised potential energy = %.9f kcal/mol"
                % final_state.getPotentialEnergy().value_in_unit(kilocalories_per_mole))
            if 1.3 <= delta_G:
                delta_g_list.append(delta_G)
            logging.info("Folding free energy  = %.9f kcal/mol"
                % delta_G)
            PDBFile.writeFile(
            simulationII.topology, final_state.getPositions(), open(pdbfolded, "w"))
    min_G = min(delta_g_list)
    ordered_G = sorted(delta_g_list)
    lowest_3 = np.average(ordered_G[:3])
    lowest_5 = np.average(ordered_G[:5])
    lowest_10 = np.average(ordered_G[:10])
    lowest_20 = np.average(ordered_G[:20])
    logging.info("Mean folding free energy (lowest 20)  = %.9f kcal/mol"
            % lowest_20)
    logging.info("Mean folding free energy (lowest 10)  = %.9f kcal/mol"
            % lowest_10)
    logging.info("Mean folding free energy (lowest 5)  = %.9f kcal/mol"
            % lowest_5)
    logging.info("Mean folding free energy (lowest 3)  = %.9f kcal/mol"
            % lowest_3)
    logging.info("Minimum folding free energy  = %.9f kcal/mol"
            % min_G)
    new_row = {'name': [stem], 'min': [min_G], 'lowest3': [lowest_3], 'lowest5': [lowest_5], 'lowest10': [lowest_10], 'lowest20': [lowest_20]}

    df_nr = pd.DataFrame(new_row)
    df_nr.to_csv('data.csv', mode='a', index=False)

    
if __name__ == '__main__':
    arguments = docopt(__doc__, version='openmm-minimise.py')
    logging.getLogger().setLevel(logging.INFO)
    main()
