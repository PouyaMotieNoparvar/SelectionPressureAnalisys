import re
import os

# Define input and output directories
input_dir = "06_tree/trees_cleaned"  # Directory with tree files to be processed
output_dir = "06_tree/trees_final"   # Directory to save fully cleaned tree files
os.makedirs(output_dir, exist_ok=True)

# Loop over all files with .tree extension in the input directory
for filename in os.listdir(input_dir):
    if filename.endswith(".tree"):
        input_path = os.path.join(input_dir, filename)
        output_path = os.path.join(output_dir, filename)

        # Read the contents of the tree file
        with open(input_path, 'r') as file:
            tree = file.read()

        # Replace very small branch lengths with 0.0
        cleaned_tree = re.sub(r':0\.000001[0-9]*', ':0.0', tree)

        # Write the cleaned tree to the output directory
        with open(output_path, 'w') as file:
            file.write(cleaned_tree)

        print(f"Cleaned tree saved as {output_path}")

print("All tree files have been processed and saved in", output_dir)
