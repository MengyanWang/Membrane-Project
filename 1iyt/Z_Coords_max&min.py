filename = 'z_coords.txt'

with open(filename, "r") as file:
    z_file = file.read()

z_value = z_file.split()

z = [float(num) for num in z_value ]

max_z = max(z)
min_z = min(z)

print(max_z, min_z)