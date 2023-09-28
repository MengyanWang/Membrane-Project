import sys
import subprocess

# get the arguments
pdb_filename = sys.argv[1]
z_change = sys.argv[2]  # a string like "+15.0" or "-15.0"
output_filename = sys.argv[3]  # the output PDB file name

# Validate z_change
if not z_change.startswith(("+", "-")) or not z_change[1:].replace('.', '', 1).isdigit():
    print("Invalid z_change argument. It should start with '+' or '-' followed by a number. python input.pdb z_change output.pdb")
    sys.exit(1)

# The TCL script
tcl_script = f"""
# Load the pdb file
mol new {pdb_filename}

# Select all atoms from the structure
set sel [atomselect 0 all]

# Get the z coordinate from the selection
set z_coords [$sel get z]

# Create empty list for new z coordinates
set new_z []

# Apply the operation to every z coordinate and append them to the array
foreach z $z_coords {{
    lappend new_z [expr $z {z_change}]
}}

# Set the new coordinates to the original structure
$sel set z $new_z

# Write the modified structure to a file
animate write pdb {output_filename} sel

# Quit VMD
quit
"""

# Write the TCL script to a file
with open('vmd_script.tcl', 'w') as f:
    f.write(tcl_script)

# Run VMD with the script file as an argument
subprocess.run(['vmd', '-dispdev', 'text', '-e', 'vmd_script.tcl'])
