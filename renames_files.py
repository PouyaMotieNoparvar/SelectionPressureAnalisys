import os

# Set the folder path where files are located
folder_path = "./02_orthogroups_cds"  # Replace with your folder path

# Loop through each file in the folder
for filename in os.listdir(folder_path):
    # Check if the filename contains "_ortho_cds"
    if "_ortho_cds" in filename:
        # Create new filename by replacing "_ortho_cds" with an empty string
        new_filename = filename.replace("_ortho_cds", "")
        
        # Define the full old and new paths
        old_path = os.path.join(folder_path, filename)
        new_path = os.path.join(folder_path, new_filename)
        
        # Rename the file
        os.rename(old_path, new_path)
        print(f'Renamed: "{filename}" to "{new_filename}"')

print("Renaming complete.")
