import sys
import mdtraj as md
input_pdb = sys.argv[1]
output_pdb = sys.argv[2]
stride = int(sys.argv[3])

trajactory = md.load(input_pdb)
subsampled_trajactory = trajactory[::stride]
subsampled_trajactory.save(output_pdb)