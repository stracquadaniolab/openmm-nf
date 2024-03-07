#!/usr/bin/env python3
"""openmm-md

Minimise the potential energy of the wildtype protein structure and mutated variants, according to AMBER14 potentials 
Usage:
openmm-md.py [--i=<pdb>] [--solv-mol=<mol>] [--temp==<temp>] [--pres=<pres>] [--nvt=<nvt>] [--npt=<npt] [--rep=<rep>]
Options:
--i=<pdb>          Input the fixed PDB files from the Mutant Maker process
--solv-mol=<mol>   Number of solvent molecules (default: solvation with H2O) 
--temp=<temp>      Temperature of MD thermostat 
--pres=<pres>      Pressure of MD barostat (NPT MD only)
--nvt=<nvt>        Number of NVT trajectory steps to initiate MD (each step = 1fs)
--npt=<npt>        Number of NPT steps in principal MD trajectory (each step = 1fs)
--rep=<rep>        Reporting rate wherefrom PDB coordinates and therodynamic conditions are outputted to trajectory (PDB) and log (CSV) files respectively
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

def setup_modeller(pdb):
    modeller = Modeller(pdb.topology, pdb.positions)
    return modeller

def setup_system(modeller, forcefield, solvmol: str):
    mol=int(solvmol)
    modeller.addSolvent(forcefield, numAdded=mol)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
    return system

def setup_simulation(modeller, system, temp: str):
    temp_fl = int(temp)
    integrator = LangevinMiddleIntegrator(temp_fl*kelvin, 1/picosecond, 0.001*picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    return simulation

def main():
    arguments = docopt(__doc__, version='openmm-md.py')
    pdb = get_pdb(arguments['--i'])
    forcefield = setup_forcefield()
    modeller = setup_modeller(pdb)
    system = setup_system(modeller, forcefield, arguments['--solv-mol'])
    simulation = setup_simulation(modeller, system, arguments['--temp'])
    print("Minimizing energy")
    simulation.minimizeEnergy()
    print("Running NVT")
    reprate = int(arguments['--rep'])
    nvt = int(arguments['--nvt'])
    npt = int(arguments['--npt'])
    pres = int(arguments['--pres'])
    temp = int(arguments['--temp'])
    str_name = arguments['--i']
    stem = str_name.replace("_opt.pdb", "")
    simulation.reporters.append(PDBReporter(stem + "_MD.pdb", reprate))
    simulation.reporters.append(StateDataReporter(stdout, reprate, step=True, potentialEnergy=True, temperature=True, volume=True))
    simulation.reporters.append(StateDataReporter(stem + "_log.csv", reprate, step=True,
        potentialEnergy=True, temperature=True, volume=True))
    simulation.step(nvt)
    system.addForce(MonteCarloBarostat(pres*bar, temp*kelvin))
    simulation.context.reinitialize(preserveState=True)
    print("Running NPT")
    simulation.step(npt)
if __name__ == '__main__':
    main()
