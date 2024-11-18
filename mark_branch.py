import os
from ete3 import Tree

def mark_ancestral_branch(newick_file, species_name="Arabidopsis_thaliana", output_folder="./06_tree/trees_xlabeled"):
    """
    Reads a Newick file, marks the ancestral branch to a given species with #1, 
    and writes the modified tree to a new file in the output folder.
    """
    # Load the Newick tree
    tree = Tree(newick_file, format=1)
    
    # Find the specified species node
    target_node = tree.search_nodes(name=species_name)
    if not target_node:
        print(f"Species '{species_name}' not found in {newick_file}. Skipping.")
        return
    target_node = target_node[0]  # Get the node itself

    # Get the ancestor node of the target species
    ancestor = target_node.up
    if ancestor:
        ancestor.add_feature("label", "#1")  # Mark the branch with #1

    # Ensure output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Save the modified tree
    output_path = os.path.join(output_folder, os.path.basename(newick_file))
    tree.write(outfile=output_path, format=1)
    print(f"Marked tree saved to {output_path}")

def process_all_trees(input_folder, species_name="Arabidopsis_thaliana"):
    """
    Processes all Newick files in a given folder, marking the ancestral branch for each tree.
    """
    for filename in os.listdir(input_folder):
        if filename.endswith(".nwk") or filename.endswith(".tree") or filename.endswith(".txt"):
            newick_file = os.path.join(input_folder, filename)
            mark_ancestral_branch(newick_file, species_name)

# Specify the folder containing Newick files
input_folder = "./06_tree/trees_plain"
process_all_trees(input_folder, species_name="Arabidopsis_thaliana")
