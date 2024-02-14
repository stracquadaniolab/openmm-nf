from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
from sys import stdout

pdb = PDBFile(sys.argv[1])
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
modeller = Modeller(pdb.topology, pdb.positions)
modeller.deleteWater()
residues=modeller.addHydrogens(forcefield)
modeller.addSolvent(forcefield, numAdded=sys.argv[4] padding=1.0*nanometer)
system = forcefield.createSystem(modeller.topology, nonbondedMethod=PME, nonbondedCutoff=1.0*nanometer, constraints=HBonds)
integrator = LangevinMiddleIntegrator(sys.argv[2]*kelvin, 1/picosecond, 0.004*picoseconds)
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)
print("Minimizing energy")
simulation.minimizeEnergy()
print("Running NVT")
simulation.reporters.append(PDBReporter(sys.argv[8], sys.argv[7]))
simulation.reporters.append(StateDataReporter(stdout, sys.argv[7], step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.reporters.append(StateDataReporter(sys.argv[9], sys.argv[7], step=True,
        potentialEnergy=True, temperature=True, volume=True))
simulation.step(sys.argv[5])
system.addForce(MonteCarloBarostat(sys.argv[3]*bar, sys.argv[2]*kelvin))
simulation.context.reinitialize(preserveState=True)
print("Running NPT")
simulation.step(sys.argv[6])
