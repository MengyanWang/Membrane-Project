import prody as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(
    description="This is a python3 script to automatic copy the template file, run simulations")
parser.add_argument("pdb_file")
args = parser.parse_args()

pdb = pd.parsePDB(args.pdb_file)
coords = pdb.select('protein and resid 1 to 17').getCoordsets()
coords = np.reshape(coords, (coords.shape[0], -1, 3))

com_list = []
com = pd.calcCenter(coords)
print(com)

with open('Center_Of_Mass_for_N17', 'a') as f:  # Use 'a' mode for append
    f.write(args.pdb_file + ':' + str(com) + '\n')  # Add '\n' at the end to insert a new line
