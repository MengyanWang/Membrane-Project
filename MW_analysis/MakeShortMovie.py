import MDAnalysis as mda
from MDAnalysis.coordinates.PDB import PDBWriter
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("step")
args = parser.parse_args()

# Replace these with the input and output file paths
input_pdb_file = "movie.pdb"
output_pdb_file = "short_movie.pdb"

# Load the input PDB file
u = mda.Universe(input_pdb_file)

# Set the frame step size
step = int(args.step)

# Create a PDB writer
with PDBWriter(output_pdb_file, multiframe=True) as pdb_writer:
    # Iterate through the trajectory, extracting every 100,000th frame
    for ts in u.trajectory[::step]:
        pdb_writer.write(u.atoms)

print("Extraction complete.")
