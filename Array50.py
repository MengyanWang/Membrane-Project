from Bio.PDB import *
from itertools import product

# Load the PDB file
parser = PDBParser()
structure = parser.get_structure('surf', 'surf.pdb')

# Assume that we have only one model in the PDB file
model = structure[0]

# Get the chain A
chain_A = model['A']

# Define the translation distance in angstroms
distance = 26

# Create a list of unique chain IDs
# Starting from 'B' as 'A' is already in the structure
chain_ids = [chr(i) for i in range(ord('B'), ord('Z')+1)] + [chr(i) for i in range(ord('a'), ord('z')+1)]
chain_ids = chain_ids[:50]

# Create a 5*10 array
for idx, (i, j) in enumerate(product(range(5), range(10))):
    # Skip the first iteration as chain A already exists
    if i == 0 and j == 0:
        continue

    # Copy chain A to a new chain
    new_chain = chain_A.copy()

    # Get a new chain id from the list
    new_chain.id = chain_ids[idx]

    # Define the translation vector
    vector = Vector(i * distance, j * distance, 0)

    # Apply the translation to every atom in the new chain
    for atom in new_chain.get_atoms():
        current_coord = atom.get_coord()
        new_coord = Vector(current_coord[0], current_coord[1], current_coord[2]) + vector
        atom.set_coord(new_coord)
        
    # Add the new chain to the model
    model.add(new_chain)

# Save the new structure to a new PDB file
io = PDBIO()
io.set_structure(structure)
io.save('50.pdb')
