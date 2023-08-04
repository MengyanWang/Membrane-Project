import os
import MDAnalysis as mda

# Define the paths of the input pdb files
pdb_paths = [f"../frag3rd/rep{i}/movie.pdb" for i in range(1, 21)]

# Define the path of the output directory
output_dir = "lastframe"

# Loop through the input pdb files and extract the last frame
for i, pdb_path in enumerate(pdb_paths):
    # Load the pdb file into an MDAnalysis Universe object
    u = mda.Universe(pdb_path)
    
    # Get the number of frames in the trajectory
    n_frames = len(u.trajectory)
    
    # Get the last frame
    last_frame = u.trajectory[-1]
    
    # Define the output file name
    output_file = f"end-{i+1}.pdb"
    
    # Define the output file path
    output_path = os.path.join(output_dir, output_file)
    
    # Write the last frame to the output file
    with mda.Writer(output_path, bonds=None, n_atoms=u.atoms.n_atoms, format="PDB") as w:
        w.write(u)
