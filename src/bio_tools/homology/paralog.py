#%%
from Bio import Phylo
from io import StringIO
from pycirclize import Circos
from pathlib import Path
from typing import Callable


from bio_tools.viz.tree import maximal_monophyletic_clades_with_singletons
from bio_tools.phylo.ete3_utils import is_fork_node
#%% ---------HELPER FUNCTIONS---------------





TEST_TREE_NEWICK = """
(   
    (
        (
            (
                (
                    (
                        (
                            F1_E:0.7
                        ):0.45,
                        (
                            (
                                F1_D:0.5
                            ):0.45,
                            (
                                (
                                    F1_C:0.6
                                ):0.7,
                                (
                                    F1_B:0.9,
                                    F1_A:0.7
                                ):0.3
                            ):0.3
                        ):0.3
                    ):0.8,
                    (
                        (
                            F2_E
                        ):0.145,
                        (
                            (
                                F2.1_C
                            ):0.145,
                            (
                                (
                                    F2.2_C
                                ):0.17,
                                (
                                    F2_B:0.09,
                                    F2_A:0.07
                                ):0.13
                            ):0.09
                        ):0.09
                    ):0.08
                ):0.4,
                (
                    F3_A:0.05,
                    F3_C:0.06
                ):0.4
            ):0.4
        ):0.4,
        (
            (
                (
                    F4_E:0.7
                ):0.45,
                (
                    (
                        F4_D:0.5
                    ):0.45,
                    (
                        (
                            F4.1_C:0.04,
                            F4.2_C:0.05
                        ):0.03,
                        (
                            F4_B:0.9,
                            F4_A:0.7
                        ):0.3
                    ):0.3
                ):0.3
            ):0.4,
            (
                (
                    (
                        (                
                            (
                                (
                                    F5.1_A:0.04,
                                    F5.1_B:0.05
                                ):0.3,
                                (
                                    F5.2_A:0.04,
                                    F5.2_B:0.05
                                ):0.3
                            ):0.3,
                            (
                                (
                                    F5.3_B 
                                ):0.2
                            ):0.23 
                        ):0.3,
                        (
                            F5_D:0.9 
                        ):0.2
                    ):0.4
                )
            ):0.35
        ):0.3
    ):0.3, 
    (
        (
            (
                (
                    (
                        F6_A
                    ):0.3, 
                    (
                        (
                            (
                                F6_B
                            ):0.3, 
                            (
                                (
                                    F6.1_C:2.0
                                ):0.92,
                                (
                                    (
                                        (
                                            F6.2_C
                                        ):0.09, 
                                        (
                                            (
                                                (
                                                    F6.3_C:0.02
                                                )0.02,
                                                (
                                                    F6.4_C:0.03
                                                )0.03
                                            )
                                        ):0.2
                                    ):0.4
                                )   
                            ):0.2
                        ):0.3
                    )   
                ):0.2 
            )
        ):0.3, 
        (
            (
                (
                    (
                        (
                            (
                                F7_A:0.3
                            )0.3,
                            (
                                F7.1_B:0.4
                            )0.2
                        ):0.2,
                        (
                            (
                                F7_C:0.4,
                                F7.2_B:0.4
                            ):0.3, 
                            (
                                F7_D:0.4
                            ):0.3
                        ):0.3
                    ):0.3
                ):0.3
            )0:3
        )
    )
);
"""
"""
tree = Phylo.read(StringIO(TEST_TREE_NEWICK), "newick")
circos, tv = Circos.initialize_from_tree(
    tree_data=tree, 
    start=60,
    end= 300,
    r_lim=(30, 100),
)


# collapse for sure
clades_collapse = maximal_monophyletic_clades_with_singletons(
    tree=tree, 
    target_leaves=["F4.2_C", "F6.3_C", "F6.4_C"])
for clade in clades_collapse:
    tv.set_node_line_props(
        clade, 
        color="red", 
        apply_label_color=True
    )

# collapse in question
clades_in_question = maximal_monophyletic_clades_with_singletons(
    tree=tree, 
    target_leaves=[
        "F2_A", "F3_A",
        "F5.1_A", "F5.2_A","F5.1_B", "F5.2_B", 
        "F6.1_C", 
        "F7.1_B", "F7.2_B"
    ])
for clade in clades_in_question:
    tv.set_node_line_props(
        clade, 
        color="orange", 
        apply_label_color=True
    )


fig = circos.plotfig()"""

# %%

from ete3 import Tree


def collect_consecutive_leaves_by_species(
    newick: str,
    species_extractor: Callable[[str], str] = lambda name: name.split("__")[-1],
) -> list[list[Tree]]:
    """
    Pre-scan a gene tree to collect consecutive leaves belonging to the same species.

    This function traverses a phylogenetic gene tree in preorder and groups
    neighboring (consecutive in traversal order) leaves that share the same
    species identifier. Each such group represents a *potential* paralog set.

    Parameters
    ----------
    newick : str
        Newick-formatted gene tree.

    species_extractor : callable, optional
        Function that extracts species information from a leaf name.
        By default, assumes leaf labels follow:
            <gene_id__any_info__species>
        and extracts the species as the last "__"-separated field.

    Returns
    -------
    List[List[str]]
        A list of groups, each group being a list of leaf names belonging
        to the same species and appearing consecutively in traversal order.

    Notes
    -----
    - This function only identifies *potential* paralogs.
    - No tree topology or duplication/speciation logic is applied here.
    - Validation must be done downstream.
    """

    tree = Tree(newick, format=1)

    potential_groups = []
    current_group = []

    previous_species = None

    for node in tree.traverse(strategy="preorder"):
        if not node.is_leaf():
            continue

        species = species_extractor(node.name)

        if species == previous_species:
            current_group.append(node)
        else:
            if len(current_group) > 1:
                potential_groups.append(current_group)

            current_group = [node]
            previous_species = species

    # Handle last group
    if len(current_group) > 1:
        potential_groups.append(current_group)

    return potential_groups


species_extractor_test_tree = lambda name: (name.split("_")[-1])
potential_groups = collect_consecutive_leaves_by_species(newick=TEST_TREE_NEWICK, species_extractor=species_extractor_test_tree)
print(potential_groups)


#%%








def _is_redundant_pair(node1, node2, max_forks=2):
    """
    Two leaves are redundant if their LCA is within max_forks
    internal forking nodes from BOTH leaves.
    """
    lca = node1.get_common_ancestor(node2)

    d1 = fork_distance(node1, lca)
    d2 = fork_distance(node2, lca)

    return d1 <= max_forks and d2 <= max_forks

def _extend_redundant_clade(leaves, max_forks=2):
    """
    Extend a redundant clade by chaining adjacency:
    each new leaf must be redundant with the previous one.
    """
    clade = [leaves[0]]

    for leaf in leaves[1:]:
        if _is_redundant_pair(leaf, clade[-1], max_forks):
            clade.append(leaf)
        else:
            break

    return clade


def detect_redundant_paralog_clades(neighbor_groups, max_forks=2):
    """
    neighbor_groups: list of lists of neighboring leaves
                     (same species, traversal-adjacent)

    Returns list of redundant clades.
    """
    redundant = []

    for group in neighbor_groups:
        i = 0
        while i < len(group) - 1:
            clade = _extend_redundant_clade(group[i:], max_forks)

            if len(clade) > 1:
                redundant.append(clade)
                i += len(clade)
            else:
                i += 1

    return redundant

# %%
