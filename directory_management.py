"""
import os
import shutil

# Define source and destination directories
source_dir = "./05_fasta2phylip/phylip"
dest_base_dir = "./07_codeml/model_0"

# Ensure the destination directory exists
os.makedirs(dest_base_dir, exist_ok=True)

# Iterate through all `.phy` files in the source directory
for file_name in os.listdir(source_dir):
    if file_name.endswith(".phy"):
        # Extract the basename (without extension)
        base_name = os.path.splitext(file_name)[0]

        # Create a corresponding folder in the destination directory
        dest_folder = os.path.join(dest_base_dir, base_name)
        os.makedirs(dest_folder, exist_ok=True)

        # Copy the file to the new folder
        src_file = os.path.join(source_dir, file_name)
        dest_file = os.path.join(dest_folder, file_name)
        shutil.copy(src_file, dest_file)

        print(f"Copied {file_name} to {dest_folder}")
"""
import os
import shutil

# Paths to the source and destination directories
source_dir = "./06_tree/bestTree"
destination_root = "./07_codeml/model_0"

# Iterate over all files in the source directory
for file_name in os.listdir(source_dir):
    if file_name.startswith("RAxML_bestTree.") and file_name.endswith(".tree"):
        # Extract the basename
        basename = file_name[len("RAxML_bestTree."):-len(".tree")]
        
        # Determine the destination folder
        destination_folder = os.path.join(destination_root, basename)
        
        if os.path.isdir(destination_folder):
            # Define the source file path and the new destination file path
            source_file = os.path.join(source_dir, file_name)
            destination_file = os.path.join(destination_folder, f"{basename}.tree")
            
            # Copy and rename the file
            shutil.copy(source_file, destination_file)
            print(f"Copied and renamed {source_file} to {destination_file}")
        else:
            print(f"Destination folder {destination_folder} does not exist, skipping.")

print("Operation completed.")
