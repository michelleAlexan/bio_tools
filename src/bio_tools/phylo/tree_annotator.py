#%%
import yaml
from ete3 import Tree, TreeStyle


PATH_CHAR_2ODD_TREE = "/Users/michellealexander/projects/bio_tools/src/bio_tools/phylo/2ODD_baits_tree"
PATH_CHAR_2ODD_ANNO_CONFIG = "/Users/michellealexander/projects/bio_tools/src/bio_tools/phylo/2ODD_annotator_config.yml"
with open(PATH_CHAR_2ODD_ANNO_CONFIG) as stream:
    y = yaml.safe_load(stream)
    colors_2ODD_funcs = y["colors_2ODD_functions"]

#%%
t = Tree(PATH_CHAR_2ODD_TREE)
ts = TreeStyle()
t.render("tree.png", tree_style=ts)

#%%
# color ranges annotator files
file_annotator_itol_func_colors = "/Users/michellealexander/projects/bio_tools/src/bio_tools/baits/itol_func_colors.txt"
with open(file_annotator_itol_func_colors, "w") as file:
    file.write("TREE_COLORS")
    file.write("SEPARATOR COMMA")
    for k, v in colors_2ODD_funcs.items():
        file.write(f"{k},{v}")


# create annotator file for characterized bait sequences
file_annotator_itol_char_2ODDs_star = "/Users/michellealexander/projects/bio_tools/src/bio_tools/baits/itol_func_colors.txt"
with open(file_annotator_itol_func_colors, "w") as file:
    file.write("LABEL")
    file.write("SEPARATOR COMMA")
    for seq in baits.ingroup_df["fasta_id"]:
        file.write(f"{seq},star")