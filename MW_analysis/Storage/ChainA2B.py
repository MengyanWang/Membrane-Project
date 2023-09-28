import sys

def change_chain(input_file, output_file):
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                chain_id = line[21]
                if chain_id == "A":
                    # Replace chain A with chain B
                    modified_line = line[:21] + "B" + line[22:]
                    outfile.write(modified_line)
                else:
                    outfile.write(line)
            else:
                outfile.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python change_chain.py <input_pdb_file> <output_pdb_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    change_chain(input_file, output_file)
