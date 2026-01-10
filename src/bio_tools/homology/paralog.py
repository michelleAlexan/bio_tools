#%%
from Bio import Phylo
from io import StringIO
from pycirclize import Circos
from pathlib import Path
from typing import Callable
from ete3 import Tree


from bio_tools.viz.tree import maximal_monophyletic_clades_with_singletons
from bio_tools.phylo.ete3_utils import is_fork_node

#%% --------- CODE TO VISUALIZE CONCEPT OF DETECTING NEIGHBORING PARALOGS ---------------
# For better understanding the concept of neighboring paralog detection, 
# plot an example tree to visualize which leaves should be detected and which
# leaves are edge cases that should NOT be detected
VISUALIZE_EXAMPLE = False

if VISUALIZE_EXAMPLE:
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
                                        F1.1_E:0.9,
                                        F1.2_E:0.7
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
                                F4.1_A:0.9,
                                F4.2_A:0.7
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
                                        F5.2_A:0.05
                                    ):0.3,
                                    (
                                        F5.3_A:0.04,
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
                                    F7.2_B:0.4,
                                    F7.3_B:0.4
                                ):0.3, 
                                (
                                    F7.4_B:0.4
                                ):0.3
                            ):0.3
                        ):0.3
                    ):0.3
                )0:3
            )
        )
    );
    """

    tree = Phylo.read(StringIO(TEST_TREE_NEWICK), "newick")
    circos, tv = Circos.initialize_from_tree(
        tree_data=tree, 
        start=60,
        end= 300,
        r_lim=(30, 100),
    )


    # neighboring paralogs to detect => green
    clades_collapse = maximal_monophyletic_clades_with_singletons(
        tree=tree, 
        target_leaves=[
            "F4.1_C", "F4.2_C",
            "F4.1_A", "F4.2_A",
            "F5.1_A", "F5.2_A",
            "F6.1_C", "F6.2_C", "F6.3_C", "F6.4_C",
            "F7.2_B", "F7.3_B", "F7.4_B",
            "F1.1_E", "F1.2_E"
            ]
        )
    for clade in clades_collapse:
        tv.set_node_line_props(
            clade, 
            color="green", 
            apply_label_color=True
        )

    # neighboring leaves of the same organism but that do not share a LCA, 
    # aka non neighboring paralogs which the algo should not return => red
    clades_in_question = maximal_monophyletic_clades_with_singletons(
        tree=tree, 
        target_leaves=[
            "F2.1_C", "F2.2_C",
            "F2_A", "F3_A",
            "F5.3_A",
            "F5.2_B", "F5.3_B",
            "F7.1_B", 
            "F2_E"
        ])
    for clade in clades_in_question:
        tv.set_node_line_props(
            clade, 
            color="red", 
            apply_label_color=True
        )


    fig = circos.plotfig()

# %% ------------ NEIGHBORING PARALOGS DETECTION -----------------


def collect_neighboring_leaves_by_species(
    tree: Tree,
    species_extractor: Callable[[str], str] = lambda name: name.split("__")[-1],
) -> list[tuple[list[Tree], Tree]]:
    """
    Pre-scan a gene tree to collect neighboring leaves belonging to the same species.

    This function traverses a phylogenetic gene tree in preorder and groups
    neighboring (consecutive in traversal order) leaves that share the same
    species identifier. Each such group represents a *potential* paralog set.

    Parameters
    ----------
    tree : ete3.Tree
        The input argument must be an ete3.Tree that has already been instantiated. 
        This is important for tracking internal node ids that are definded by random integers. 

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
    - No Last Common Ancestor logic is applied here.
    - Validation must be done downstream.
    """


    results: list[tuple[list[Tree], Tree]] = []

    current_group: list[Tree] = []
    current_lca: Tree | None = None
    previous_species = None

    for node in tree.traverse(strategy="preorder"):
        if not node.is_leaf():
            continue

        species = species_extractor(node.name)

        if species == previous_species:
            # extend current group
            current_group.append(node)

            # update running LCA
            if current_lca is None:
                # group has exactly 2 leaves now
                current_lca = current_group[0].get_common_ancestor(node)
            else:
                current_lca = current_lca.get_common_ancestor(node)

        else:
            # finalize previous group
            if len(current_group) > 1:
                results.append((current_group, current_lca))

            # start new group
            current_group = [node]
            current_lca = None
            previous_species = species

    # handle last group
    if len(current_group) > 1:
        results.append((current_group, current_lca))

    return results


def clade_is_species_homogeneous(lca: Tree, species: str, 
                                 species_extractor: Callable[[str], str] = lambda name: name.split("__")[-1],) -> bool:
    """
    Return True if all leaves under lca belong to the same species.
    """
    for leaf in lca.iter_leaves():
        if species_extractor(leaf.name) != species:
            return False
    return True


# %%
def true_redundant_paralog_clades(
        groups_of_neighboring_leaves_of_same_species: list[tuple[list[Tree], Tree]],
        species_extractor: Callable[[str], str] = lambda name: name.split("__")[-1],
        
    ) -> list[tuple[list[Tree], Tree]]:

    results: list[tuple[list[Tree], Tree]] = []

    for group, lca in groups_of_neighboring_leaves_of_same_species:
        species = species_extractor(group[0].name)

        # if the given group contains only leaves, 
        # that are redundant paralogs, append them to the result list
        lca_clade_is_species_homogeneous = clade_is_species_homogeneous(lca, species, species_extractor)
        if lca_clade_is_species_homogeneous:
            results.append((group, lca))

        # otherwise: oh boy, here we go:
        else:
            # if the group only consists of 2 members and the clade under the lca-node of the two neighboring leaves (remember: same species) 
            # is not species-homogenous, then discard the group
            if len(group) == 2: 
                continue

            # there may be multiple clades within a group that are redundant paralogs
            # this is why we need a list (group_results) to keep track
            group_results: list[tuple[list[Tree], Tree]] = []
            done = False

            # Find all clades that are redundant paralogs within the group
            # expand each clade as long the invariant:
            #  "all leaves under the current lca belong to the same species" holds true
            # for this, iterate over the leaf nodes within the group
            # start with the first node...
            for i in range(len(group) - 1):
                # ... and compare to all other leaves one by one
                true_members: list[Tree] = [group[i]]
                for j in range(i + 1, len(group)):

                    # Retrieve the lca to then check if the invariant still holds true.
                    # If so, move on to the next leaf.
                    # If the current leaf (group[j]) is already the last leaf of the group, 
                    # handle everthing else in the else condition. 
                    lca_temp = group[i].get_common_ancestor(group[j])
                    if (clade_is_species_homogeneous(lca_temp, species, species_extractor) 
                        and j < len(group) - 1):
                        true_members.append(group[j])
                        continue
                    
                    # here, we reach the end of the (sub-) group of redundant paralogs
                    # there are differnet cases to be aware of
                    else:
                        # The current range of leaves must span more than one leaf 
                        # - otherwise it wouldnt be a group.
                        # In example test tree, this exludes the F7.1_B single group
                        if j - i > 1:

                            # Example: there are 5 members in the group. 
                            # Members at index 0 to 2 belong to one true redundant paralogous clade 
                            # (invariant holds true up to index 2).
                            # This is why we are now in the else condition with j being 3. 
                            # If there are at least 2 more members in the group, 
                            # we have to proceed later with updated indices i and j
                            if j < len(group) - 2:
                                # not the current, but the previous j (j-1) closes the true redundant paralogous clade 
                                # get the lca of the group members and collect tree nodes
                                lca_temp = group[i].get_common_ancestor(group[j - 1]) 
                                if not clade_is_species_homogeneous(lca_temp, species, species_extractor):
                                    raise ValueError("You are creating a 'true redundant paralogous clade', " \
                                    "but at the same time you are violating the invariant that all leaves under the lca must come from the same species only")
                                
                                redundant_paralogous_clade = (true_members, lca_temp)
                                group_results.append(redundant_paralogous_clade)

                                # one true redundant paralogous subgroup within the group is closed
                                # since there are at least 2 more tree nodes we have to check, update i and j
                                i = j
                                j = j + 1
                                break

                            # the last member of the group has been reached. 
                            else: 
                                # there are two cases now:
                                # first: the last member of the group does NOT belong to redundant_paralogous_clade
                                lca_temp = group[i].get_common_ancestor(group[j])
                                if not clade_is_species_homogeneous(lca_temp, species, species_extractor):
                                    j_last_true_member = j - 1

                                # second: the last member DOES belong to the redundant_paralogous_clade
                                else:
                                    j_last_true_member = j 
                                    true_members.append(group[j])
                                    
                                lca_temp = group[i].get_common_ancestor(group[j_last_true_member]) 
                                if not clade_is_species_homogeneous(lca_temp, species, species_extractor):
                                    raise ValueError("You are creating a 'true redundant paralogous clade', " \
                                    "but at the same time you are violating the invariant that all leaves under the lca must come from the same species only")
                                
                                redundant_paralogous_clade = (true_members, lca_temp)
                                group_results.append(redundant_paralogous_clade)
                        else:
                            break

                if done:
                    break  
            if group_results:        
                results.extend(group_results)
    for group, lca in results:
        print(lca)

    return results











                                

#%%








# %%
