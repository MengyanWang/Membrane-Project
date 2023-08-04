
# Load the pdb file
mol new input.pdb

# Select all atoms from the structure
set sel [atomselect 0 all]

# Get the z coordinate from the selection
set z_coords [$sel get z]

# Create empty list for new z coordinates
set new_z []

# Apply the operation to every z coordinate and append them to the array
foreach z $z_coords {
    lappend new_z [expr $z +14.8400566]
}

# Set the new coordinates to the original structure
$sel set z $new_z

# Write the modified structure to a file
animate write pdb surf.pdb sel

# Quit VMD
quit
