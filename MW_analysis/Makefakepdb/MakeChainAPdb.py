#!/usr/bin/env python3

import os
import sys

squence_file = sys.argv[1] + '.seq'
pml_file = sys.argv[2] + '.pml'
pdb_file = sys.argv[3] + '.pdb'

with open(squence_file,'r') as fopen:
    lines = fopen.readlines()
    for line1 in lines:
       line = line1.rstrip('\n')
       break 

with open (pml_file,'w') as fwrite:
    data = ''
    data += 'sequence = \"' + line.lower() + '\" \n'
    data += 'for aa in sequence: cmd._alt(aa)\n'
    data += 'alter (all),resi=str(int(resi)-1)\n'
    data += 'alter (all), chain="A"\n' # Set chain to A for all atoms
#    data += 'ray 1000 1000\n' 
    data += 'save ' + pdb_file + '\n'
    data += 'quit \n'
    fwrite.writelines(data)

os.system("pymol %s"%(pml_file))
