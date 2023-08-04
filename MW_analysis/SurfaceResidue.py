from VectorAlgebra import *
from Bio.PDB import PDBParser

def checkIfNative(xyz_CAi, xyz_CAj):
    v = vector(xyz_CAi, xyz_CAj)
    r = vabs(v)
    if r< 9.5: return True
    else: return False

# Load the structure
p = PDBParser()
structure = p.get_structure('protein', 'standardFrames//end-1_s.pdb')

# List of residues
residues = list(structure.get_residues())

# Initialize a dictionary to hold the count of contacts for each residue
contacts = {res: 0 for res in residues}

# Loop over every pair of residues
for i in range(len(residues)):
    for j in range(i+1, len(residues)):
        res_i = residues[i]
        res_j = residues[j]
        
        # Get the CA atoms
        ca_i = res_i['CA'].get_coord()
        ca_j = res_j['CA'].get_coord()
        
        # Check if the residues are in contact
        if checkIfNative(ca_i, ca_j):
            contacts[res_i] += 1
            contacts[res_j] += 1

# Filter out the surface residues (those with less than 20 contacts)
surface_residues = [res for res, count in contacts.items() if count < 20]

# Print out the surface residues
print('Surface Residues:')
for res in surface_residues:
    print(f"{res.get_id()[1]}",end=' ')
#     print(f"Surface residue: {res.get_resname()} {res.get_id()[1]}")

# Define sets of charged residues
positive_residues = {'ARG', 'LYS', 'HIS'}
negative_residues = {'ASP', 'GLU'}

# Find the intersection of all sets of surface residues
# surface_residues = set.intersection(*surface_residues_sets)

# Initialize lists for positive and negative residues
positive_surface_residues = []
negative_surface_residues = []

# Classify common surface residues as positive or negative
for res in surface_residues:
    if res.get_resname() in positive_residues:
        positive_surface_residues.append(res)
    elif res.get_resname() in negative_residues:
        negative_surface_residues.append(res)

# Print out the common positive and negative surface residues
print("Positively charged surface residues:")
for res in positive_surface_residues:
    print(res.get_id()[1],end=" ")

print("\nNegatively charged surface residues:")
for res in negative_surface_residues:
    print(res.get_id()[1],end=" ")