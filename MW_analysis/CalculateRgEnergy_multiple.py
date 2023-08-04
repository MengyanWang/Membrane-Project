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
parser.add_argument("-f", "--forces", default="forces_setup.py")
args = parser.parse_args()

original_path = os.getcwd()
setupFolderPath = os.path.dirname(args.forces)

forceSetupFile = None if args.forces is None else os.path.abspath(args.forces)

# proteinName = pdb_id = os.path.basename(args.protein)

with open("Rg&Energy.dat", "w") as file:
    file.write("Step Rg Total\n")

for i in range(1,21):

    pdb_file = f"standardFrames//end-{i}_s.pdb"
    pdbPath = os.path.abspath(pdb_file)
    pdb = md.load(pdbPath,stride=1)
    os.chdir(setupFolderPath)

    input_pdb_filename = "ews_extended-openmmawsem.pdb"

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
    for step in range(len(pdb)):
        # print(step)
        simulation.context.setPositions(pdb.openmm_positions(step))
        
    # forceGroupTable = {"Rg":2,"Total":list(range(11, 32))}
    forceGroupTable = {"Backbone":20, "Rama":21, "Contact":22, "Fragment":23, "Membrane":24, "ER":25, "TBM_Q":26, "Beta":27, "Pap":28, "Helical":29,
            "Q":1, "Rg":2, "Qc":3,
                    "Helix_orientation":18, "Pulling":19,
                    "Total":list(range(11, 32))
                    # , "Q_wat":4, "Q_mem":5, "Debye_huckel":30
                   }

    # g = set(forceGroupTable[term])

    state = simulation.context.getState(getEnergy=True, groups={forceGroupTable["Rg"]})
    rg = state.getPotentialEnergy().value_in_unit(kilojoule_per_mole)

    state = simulation.context.getState(getEnergy=True, groups=set(forceGroupTable["Total"]))
    total = state.getPotentialEnergy().value_in_unit(kilocalories_per_mole)

    # print("Rg:", rg)

    with open("Rg&Energy.dat", "a") as file:
        file.write(str(i) + " " + str(rg)+ " " + str(total) + "\n")

    # new_path = f"/home/mw88/research/ews_extended/fuse/2nd_Analysis"
    os.chdir(original_path)