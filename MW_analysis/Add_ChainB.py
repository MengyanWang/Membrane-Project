from Bio.PDB import *
import argparse

# Argument parser
parser = argparse.ArgumentParser(description='Process a PDB file.')
parser.add_argument('input_pdb', help='Input PDB file')
parser.add_argument('output_pdb', help='Output PDB file')
args = parser.parse_args()

# Load the PDB file
parser = PDBParser()
structure = parser.get_structure('surf', args.input_pdb)

# Assume that we have only one model in the PDB file
model = structure[0]

# Get the chain A
chain_A = model['A']

# Copy chain A to chain B
chain_B = chain_A.copy()

# Update the id to B
chain_B.id = 'B'

# Add the new chain to the model
model.add(chain_B)

# Define the translation vector for 25 angstrom along the x axis
vector = Vector(25, 0, 0)

# Apply the translation to every atom in the chain
for atom in chain_B.get_atoms():
    current_coord = atom.get_coord()
    new_coord = Vector(current_coord[0], current_coord[1], current_coord[2]) + vector
    atom.set_coord(new_coord)

# Save the new structure to a new PDB file
io = PDBIO()
io.set_structure(structure)
io.save(args.output_pdb)
