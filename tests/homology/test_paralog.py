#%%
import pytest
from ete3 import Tree

from bio_tools.homology.paralog import (
    collect_neighboring_leaves_by_species,
    clade_is_species_homogeneous,
    true_redundant_paralog_clades    
)

def node_names(result):
    return [([n.name for n in group], lca) for group, lca in result]

TEST_TREE_NEWICK_STR = """
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
t = Tree(TEST_TREE_NEWICK_STR)


lca_F1_1_E__F2_E = t.get_common_ancestor('F1.1_E', 'F2_E')
lca_F2_1_C__F2_2_C = t.get_common_ancestor('F2.1_C', 'F2.2_C')
lca_F2_A__F3_A     = t.get_common_ancestor('F2_A', 'F3_A')
lca_F4_1_C__F4_2_C = t.get_common_ancestor('F4.1_C', 'F4.2_C')
lca_F4_1_A__F5_3_A   = t.get_common_ancestor('F4.1_A', 'F5.3_A')
lca_F5_2_B__F5_3_B = t.get_common_ancestor('F5.2_B', 'F5.3_B')
lca_F6_Cs = t.get_common_ancestor(
    'F6.1_C', 'F6.2_C', 'F6.3_C', 'F6.4_C'
)
lca_F7_Bs = t.get_common_ancestor(
    'F7.1_B', 'F7.2_B', 'F7.3_B', 'F7.4_B'
)
@pytest.mark.parametrize(
        "lca, species, expected",
        [
            (
                lca_F1_1_E__F2_E, "E", False
            ),
            (
                lca_F2_1_C__F2_2_C, "C", False
            ),
            (
                lca_F2_A__F3_A, "A", False
            ),
            (
                lca_F4_1_C__F4_2_C, "C", True
            ), 
            (
                lca_F4_1_A__F5_3_A, "A", False
            ),        
            (
                lca_F5_2_B__F5_3_B, "B", False
            ),    
            (
                lca_F6_Cs, "C", True
            ),    
            (
                lca_F7_Bs, "B", False
            ),    
        ]
)
def test_clade_is_species_homogeneous(lca, species, expected):
    species_extractor_test_tree = lambda name: name.split("_")[-1]
    result = clade_is_species_homogeneous(
        lca=lca, 
        species=species, 
        species_extractor=species_extractor_test_tree
        )
    assert result == expected
    

def test_collect_neighboring_leaves_by_species():


    species_extractor_test_tree = lambda name: name.split("_")[-1]

    result = collect_neighboring_leaves_by_species(
        tree=t,
        species_extractor=species_extractor_test_tree
    )

    expected = [
        (['F1.1_E', 'F1.2_E', 'F2_E'], lca_F1_1_E__F2_E),
        (['F2.1_C', 'F2.2_C'], lca_F2_1_C__F2_2_C),
        (['F2_A', 'F3_A'],     lca_F2_A__F3_A),
        (['F4.1_C', 'F4.2_C'], lca_F4_1_C__F4_2_C),
        (['F4.1_A', 'F4.2_A', 'F5.1_A', 'F5.2_A', 'F5.3_A'],   lca_F4_1_A__F5_3_A),
        (['F5.2_B', 'F5.3_B'], lca_F5_2_B__F5_3_B),
        (['F6.1_C', 'F6.2_C', 'F6.3_C', 'F6.4_C'], lca_F6_Cs),
        (['F7.1_B', 'F7.2_B', 'F7.3_B', 'F7.4_B'], lca_F7_Bs),
    ]

    assert node_names(result) == expected



lca_F1_1_E__F1_2_E = t.get_common_ancestor('F1.1_E', 'F1.2_E')
lca_F4_1_A__F4_2_A   = t.get_common_ancestor('F4.1_A', 'F4.2_A')
lca_F5_1_A__F5_2_A   = t.get_common_ancestor('F5.1_A', 'F5.2_A')
lca_F7_2_B__F7_4_B = t.get_common_ancestor(
    'F7.2_B', 'F7.3_B', 'F7.4_B'
)
def test_true_redundant_paralog_clades():
    species_extractor_test_tree = lambda name: name.split("_")[-1]

    potential_paralogous_clades = collect_neighboring_leaves_by_species(
        tree=t,
        species_extractor=species_extractor_test_tree
    )

    result = true_redundant_paralog_clades(potential_paralogous_clades, 
                                           species_extractor_test_tree)

    expected = [
        (['F1.1_E', 'F1.2_E'], lca_F1_1_E__F1_2_E),
        (['F4.1_C', 'F4.2_C'], lca_F4_1_C__F4_2_C),
        (['F4.1_A', 'F4.2_A'], lca_F4_1_A__F4_2_A),
        (['F5.1_A', 'F5.2_A'], lca_F5_1_A__F5_2_A),
        (['F6.1_C', 'F6.2_C', 'F6.3_C', 'F6.4_C'], lca_F6_Cs),
        (['F7.2_B', 'F7.3_B', 'F7.4_B'], lca_F7_2_B__F7_4_B),
    ]

    assert node_names(result) == expected


# %%
