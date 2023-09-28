import prody as pd
import numpy as np
import argparse

# pdb = pd.parsePDB("//home/mw88//research//n17 simulaiton//below.pdb")

parser = argparse.ArgumentParser(
    description="This is a python3 script to\
    automatic copy the template file, \
    run simulations")
parser.add_argument("pdb_file")
args = parser.parse_args()

pdb = pd.parsePDB(args.pdb_file)

coords = pdb.getCoordsets()

coords = np.reshape(coords, (coords.shape[0],-1,3))

com_list = []

com = pd.calcCenter(coords)

print(com)

with open('Center_Of_Mass_for_start_structure.txt', 'a') as f:
    # f.write(args.pdb_file +':',str(com))
    f.write(args.pdb_file + ':' + str(com) + '\n')
