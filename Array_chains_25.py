from Bio.PDB import *

# Load the PDB file
parser = PDBParser()
structure = parser.get_structure('surf', 'surf.pdb')

# Assume that we have only one model in the PDB file
model = structure[0]

# Get the chain A
chain_A = model['A']

# Define the translation distance in angstroms
distance = 26

# ASCII code for 'B'
chain_id = 66

# Create a 5*5 array
for i in range(5):
    for j in range(5):
        # Skip the first iteration as chain A already exists
        if i == 0 and j == 0:
            continue

        # Copy chain A to a new chain
        new_chain = chain_A.copy()

        # Generate a new chain id
        new_chain.id = chr(chain_id)

        # Define the translation vector
        vector = Vector(i * distance, j * distance, 0)

        # Apply the translation to every atom in the new chain
        for atom in new_chain.get_atoms():
            current_coord = atom.get_coord()
            new_coord = Vector(current_coord[0], current_coord[1], current_coord[2]) + vector
            atom.set_coord(new_coord)
        
        # Add the new chain to the model
        model.add(new_chain)

        # Increment the chain id
        chain_id += 1

# Save the new structure to a new PDB file
io = PDBIO()
io.set_structure(structure)
io.save('new_array.pdb')
