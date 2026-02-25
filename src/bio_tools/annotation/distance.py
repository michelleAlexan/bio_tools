#%%

from pathlib import Path
from ete4 import Tree, PhyloTree


TREE_PATH = Path("/Users/michellealexander/Documents/ingroup_outgroup/B3/post_B3_tree.nwk")
t = PhyloTree(open(TREE_PATH), sp_naming_function=lambda name: name.split('__')[-1])
tax2names, tax2lineages, tax2rank = t.annotate_ncbi_taxa(taxid_attr='species')
t.ladderize()

#%%
t.explore()
# %%
distance_matrix = t.distance_matrix()

# %%
