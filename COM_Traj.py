import prody as pd
import numpy as np

pdb = pd.parsePDB('movie.pdb')

# Get the coordinates for each frame
coords = pdb.getCoordsets()

# Reshape the array to have the shape (n_frames, n_atoms, 3)
coords = np.reshape(coords, (coords.shape[0],-1,3))

com_list = []

# Iterate over each frame in the PDB file
for frame_coords in coords:
    # compute the center of mass for the frame
    com = pd.calcCenter(frame_coords)
    # add the center of mass to the list
    com_list.append(com)
    
# print the coordinates of the center of mass for each frame
# for i,com in enumerate(com_list):
#     print(f"COM for frame {i+1}: {com}")

with open('Center_Of_Mass_for_Traj.txt', 'w') as f:
    for i,com in enumerate(com_list):
        f.write(f"COM for frame {i+1}: {com}\n")