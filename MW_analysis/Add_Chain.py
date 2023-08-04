from Bio.PDB import PDBParser, PDBIO
import argparse

# Argument parser
parser = argparse.ArgumentParser(description='Modify the chain identifier of a PDB file.')
parser.add_argument('input_pdb', help='Input PDB file')
parser.add_argument('output_pdb', help='Output PDB file')
args = parser.parse_args()

# Load the PDB file
parser = PDBParser()
structure = parser.get_structure('PDB', args.input_pdb)

# Iterate over each model in the structure
for model in structure:
    # Iterate over each chain in the model
    for chain in model:
        # Set the chain identifier to 'A'
        chain.id = 'A'

# Save the new structure to a new PDB file
io = PDBIO()
io.set_structure(structure)
io.save(args.output_pdb)
