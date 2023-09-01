import os
import argparse
import sys
import importlib.util

try:
    OPENAWSEM_LOCATION = os.environ["OPENAWSEM_LOCATION"]
    sys.path.append(OPENAWSEM_LOCATION)
except KeyError:
    print("Please set the environment variable name OPENAWSEM_LOCATION.\n Example: export OPENAWSEM_LOCATION='YOUR_OPENAWSEM_LOCATION'")
    exit()

from openmmawsem import *
from helperFunctions.myFunctions import *

parser = argparse.ArgumentParser(
    description="This code computes the radius of gyration (Rg) for a single frame PDB file."
)
parser.add_argument("protein", help="The name of the protein")
parser.add_argument("-f", "--forces", default="forces_setup.py")
args = parser.parse_args()

setupFolderPath = os.path.dirname(args.forces)
forceSetupFile = None if args.forces is None else os.path.abspath(args.forces)

proteinName = pdb_id = os.path.basename(args.protein)
pdb_file = f"{pdb_id}.pdb"
pdbPath = os.path.abspath(pdb_file)
pdb = md.load(pdbPath,stride=1)
os.chdir(setupFolderPath)


input_pdb_filename = f"{pdb_id}-openmmawsem.pdb"



chain = getAllChains("crystal_structure.pdb")
simulation_platform = "CPU"
platform = Platform.getPlatformByName(simulation_platform)

oa = OpenMMAWSEMSystem(input_pdb_filename, chains=chain, k_awsem=1.0, xml_filename=f"{OPENAWSEM_LOCATION}/awsem.xml",includeLigands=False) 

spec = importlib.util.spec_from_file_location("forces", forceSetupFile)
# print(spec)
forces = importlib.util.module_from_spec(spec)
spec.loader.exec_module(forces)
forces = forces.set_up_forces(oa, computeQ=True, submode=3, contactParameterLocation=".")
oa.addForcesWithDefaultForceGroup(forces)

collision_rate = 5.0 / picoseconds
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 2*femtoseconds)
simulation = Simulation(oa.pdb.topology, oa.system, integrator, platform)

# print("Number of atoms in PDB:", pdb.topology.n_atoms)
# print("Number of particles in the system:", simulation.system.getNumParticles())


for step in range(len(pdb)):
    print(step)
    simulation.context.setPositions(pdb.openmm_positions(step))
    
forceGroupTable = {"Rg":2}
state = simulation.context.getState(getEnergy=True, groups={forceGroupTable["Rg"]})

rg = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)
print("Rg:", rg)