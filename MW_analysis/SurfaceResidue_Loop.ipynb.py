from VectorAlgebra import *
from Bio.PDB import PDBParser
import os

def checkIfNative(xyz_CAi, xyz_CAj):
    v = vector(xyz_CAi, xyz_CAj)
    r = vabs(v)
    if r< 9.5: return True
    else: return False

def protonation_state(pH, pKa):
        return 1 / (1 + 10 ** (pH - pKa))

pH = 7.4  # set this to the appropriate pH
pKa_values = {'ARG': 12.5, 'LYS': 10.5, 'HIS': 6.0, 'ASP': 3.9, 'GLU': 4.1}

def analyze_pdb(file_name, replica_number):
    # Load the structure
    p = PDBParser()
    structure = p.get_structure('protein', file_name)

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
    # surface_residues = [res.get_id()[1] for res, count in contacts.items() if count < 20]
    surface_residues = [res for res, count in contacts.items() if count < 20]
    surface_residues_id = [res.get_id()[1] for res in surface_residues]




    # Define sets of charged residues
    positive_residues = {'ARG', 'LYS', 'HIS'}
    negative_residues = {'ASP', 'GLU'}

    # Initialize lists for positive and negative residues
    positive_surface_residues = []
    negative_surface_residues = []

    # # Classify common surface residues as positive or negative
    # for res in structure.get_residues():
    #     if res.get_resname() in positive_residues and res.get_id()[1] in surface_residues:
    #         positive_surface_residues.append(res.get_id()[1])
    #     elif res.get_resname() in negative_residues and res.get_id()[1] in surface_residues:
    #         negative_surface_residues.append(res.get_id()[1])
    # Classify common surface residues as positive or negative
    for res in surface_residues:
        if res.get_resname() in positive_residues:
            positive_surface_residues.append(res)
        elif res.get_resname() in negative_residues:
            negative_surface_residues.append(res)

    # Get the residue ids for output
    positive_residue_ids = [res.get_id()[1] for res in positive_surface_residues]
    negative_residue_ids = [res.get_id()[1] for res in negative_surface_residues]

    surface_charge = 0
    for res in positive_surface_residues:
        resname = res.get_resname()
        surface_charge += protonation_state(pH, pKa_values[resname])

    for res in negative_surface_residues:
        resname = res.get_resname()
        surface_charge -= protonation_state(pH, pKa_values[resname])

    # Define sets of hydrophobic residues
    hydrophobic_residues = {'ALA', 'ILE', 'LEU', 'MET', 'VAL', 'PHE', 'PRO', 'TRP'}
    hydrophilic_residues = {'SER', 'THR', 'ASN', 'GLN', 'ARG', 'HIS', 'LYS', 'ASP', 'GLU'}

    hydrophobic_surface_residues = []
    hydrophilic_surface_residues = []

    for res in structure.get_residues():
        if res.get_resname() in hydrophobic_residues and res.get_id()[1] in surface_residues_id:
            hydrophobic_surface_residues.append(res.get_id()[1])
        elif res.get_resname() in hydrophilic_residues and res.get_id()[1] in surface_residues_id:
            hydrophilic_surface_residues.append(res.get_id()[1])

    # Calculate the surface charge
    # surface_charge = len(positive_surface_residues) - len(negative_surface_residues)

    # Format the results
    # results = f"Replica {replica_number}\nSurface Residues: {','.join(map(str, surface_residues))}\nPositively Charged Surface Residues: {','.join(map(str, positive_surface_residues))}\nNegatively Charged Surface Residues: {','.join(map(str, negative_surface_residues))}\n\n"
    # results = f"Replica {replica_number}\nSurface Residues: {','.join(map(str, surface_residues))}\nPositively Charged Surface Residues: {','.join(map(str, positive_surface_residues))}\nNegatively Charged Surface Residues: {','.join(map(str, negative_surface_residues))}\nHydrophobic Surface Residues: {','.join(map(str, hydrophobic_surface_residues))}\nHydrophilic Surface Residues: {','.join(map(str, hydrophilic_surface_residues))}\n\n"
    results = f"Replica {replica_number}\nSurface Residues: {','.join(map(str, surface_residues_id))}\nPositively Charged Surface Residues: {','.join(map(str, positive_residue_ids))}\nNegatively Charged Surface Residues: {','.join(map(str, negative_residue_ids))}\nHydrophobic Surface Residues: {','.join(map(str, hydrophobic_surface_residues))}\nHydrophilic Surface Residues: {','.join(map(str, hydrophilic_surface_residues))}\nSurface Charge: {surface_charge}\n\n"


    # Return the formatted results
    return results

# Specify the directory containing the pdb files
directory = 'standardFrames//'

# Loop over all the PDB files
results = []
for i in range(1, 21):
    file_name = os.path.join(directory, f'end-{i}_s.pdb')
    result = analyze_pdb(file_name, i)
    results.append(result)

# Save results to a .dat file
with open('results.dat', 'w') as f:
    for result in results:
        f.write(result)
