import prody as pd
import numpy as np
import argparse
import subprocess

def write_tcl_script(pdb_filename, z_change, output_filename):
    tcl_script = f"""
    mol new {pdb_filename}
    set sel [atomselect 0 all]
    set z_coords [$sel get z]
    set new_z []
    foreach z $z_coords {{
        lappend new_z [expr $z + {z_change}]
    }}
    $sel set z $new_z
    animate write pdb {output_filename} sel
    quit
    """
    with open('vmd_script.tcl', 'w') as f:
        f.write(tcl_script)

parser = argparse.ArgumentParser(description="This is a python3 script to modify the PDB structure so that its center of mass has z=15.")
parser.add_argument("pdb_file", help="Input PDB file")
parser.add_argument("output_file", help="Output PDB file")

args = parser.parse_args()

# Parsing the PDB file
pdb = pd.parsePDB(args.pdb_file)

# Calculating center of mass
coords = pdb.getCoordsets()
coords = np.reshape(coords, (coords.shape[0],-1,3))
com = pd.calcCenter(coords)
com = np.squeeze(com)
print(f"Original Center of Mass: {com}")

# Calculating the shift required to set z component of COM to 15
z_shift = 15 - com[2]
print(f"Shift required in z-direction: {z_shift}")

# Write the TCL script to a file
write_tcl_script(args.pdb_file, z_shift, args.output_file)

# Run VMD with the TCL script
subprocess.run(['vmd', '-dispdev', 'text', '-e', 'vmd_script.tcl'])

# Parsing the new, shifted PDB file
pdb_shifted = pd.parsePDB(args.output_file)

# Calculating new center of mass
coords_shifted = pdb_shifted.getCoordsets()
coords_shifted = np.reshape(coords_shifted, (coords_shifted.shape[0],-1,3))
com_shifted = pd.calcCenter(coords_shifted)
com_shifted = np.squeeze(com_shifted)
print(f"Shifted Center of Mass: {com_shifted}")
