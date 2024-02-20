#!/usr/bin/env python3
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
from sys import stdout
from pdbfixer import PDBFixer

#pdb = PDBFile(sys.argv[1])
pdb = PDBFixer(sys.argv[1])
pdb.addMutations(mutationResidues='19', newResidues='ALA')
pdb.findMissingResidues()
pdb.findNonstandardResidues()
pdb.replaceNonstandardResidues()
pdb.removeHeterogens(True)
pdb.findMissingAtoms()
pdb.addMissingAtoms()
pdb.addMissingHydrogens(7.0)
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
#modeller.deleteWater()
#residues=modeller.addHydrogens(forcefield)
solvmol=int(sys.argv[4])
modeller.addSolvent(forcefield, numAdded=solvmol)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
temp=int(sys.argv[2])
integrator = LangevinMiddleIntegrator(temp*kelvin, 1/picosecond, 0.001*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
print("Minimizing energy")
simulation.minimizeEnergy()
print("Running NVT")
reprate=int(sys.argv[7])
simulation.reporters.append(PDBReporter(sys.argv[8], reprate))
simulation.reporters.append(StateDataReporter(stdout, reprate, step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReporter(sys.argv[9], reprate, step=True,
        potentialEnergy=True, temperature=True, volume=True))
nvtsteps=int(sys.argv[5])
pres=int(sys.argv[3])
simulation.step(nvtsteps)
system.addForce(MonteCarloBarostat(pres*bar, temp*kelvin))
simulation.context.reinitialize(preserveState=True)
print("Running NPT")
nptsteps=int(sys.argv[6])
simulation.step(nptsteps)
