#%%
import pytest
from ete4 import Tree
from pathlib import Path
import json
from io import StringIO

from bio_tools.homology.paralog import (
    collect_neighboring_leaves_by_species,
    clade_is_species_homogeneous,
    true_redundant_paralog_clades,
    detect_redundant_paralog_clades,    
    map_representative_paralog_to_all_redundant_paralogs,
    reduce_seq_collection_to_non_redundancy,
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
                                                ):0.02,
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
            ):0.3
        )
    )
);
"""
TEST_TREE_NEWICK_STR = TEST_TREE_NEWICK_STR.replace("\n", "").replace(" ", "")
t = Tree(TEST_TREE_NEWICK_STR)
species_extractor_test_tree = lambda name: name.split("_")[-1]


lca_F1_1_E__F2_E = t.common_ancestor('F1.1_E', 'F2_E')
lca_F2_1_C__F2_2_C = t.common_ancestor('F2.1_C', 'F2.2_C')
lca_F2_A__F3_A     = t.common_ancestor('F2_A', 'F3_A')
lca_F4_1_C__F4_2_C = t.common_ancestor('F4.1_C', 'F4.2_C')
lca_F4_1_A__F5_3_A   = t.common_ancestor('F4.1_A', 'F5.3_A')
lca_F5_2_B__F5_3_B = t.common_ancestor('F5.2_B', 'F5.3_B')
lca_F6_Cs = t.common_ancestor(
    'F6.1_C', 'F6.2_C', 'F6.3_C', 'F6.4_C'
)
lca_F7_Bs = t.common_ancestor(
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
    result = clade_is_species_homogeneous(
        lca=lca, 
        species=species, 
        species_extractor=species_extractor_test_tree
        )
    assert result == expected
    

def test_collect_neighboring_leaves_by_species():

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



lca_F1_1_E__F1_2_E = t.common_ancestor('F1.1_E', 'F1.2_E')
lca_F4_1_A__F4_2_A   = t.common_ancestor('F4.1_A', 'F4.2_A')
lca_F5_1_A__F5_2_A   = t.common_ancestor('F5.1_A', 'F5.2_A')
lca_F7_2_B__F7_4_B = t.common_ancestor(
    'F7.2_B', 'F7.3_B', 'F7.4_B'
)
def test_true_redundant_paralog_clades():

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
    print("result")
    for r in result:
        print(r)
    assert node_names(result) == expected


def test_detect_redundant_paralog_clades():
    expected = [
        ['F1.1_E', 'F1.2_E'],
        ['F4.1_C', 'F4.2_C'],
        ['F4.1_A', 'F4.2_A'],
        ['F5.1_A', 'F5.2_A'], 
        ['F6.1_C', 'F6.2_C', 'F6.3_C', 'F6.4_C'],
        ['F7.2_B', 'F7.3_B', 'F7.4_B'], 
    ]
    result = detect_redundant_paralog_clades(
        tree=t, 
        species_extractor=species_extractor_test_tree,
        return_as_strings=True
        )

    assert result == expected



# %%


def test_map_representative_paralog_to_all_redundant_paralogs(tmp_path):

    # mock fasta
    fasta_content = """>F4.3_C
AAAAAAAAA
>unrelated_header
QQQQQQQQQQ
>F1.1_E
EEEEEEEEEEEE
>unrelated_header2
QQQQQQ
>F1.2_E
EEEEEEEEEEEEEEEEEEE
>F4.1_C
AAAAAAA
>F4.2_C
AAAAAAAAAAAA
>unrelated_header3
QQQQQQQQ
>F5_C
AAA
"""
    example_fasta = tmp_path / "example.fasta"
    example_fasta.write_text(fasta_content)

    example_grouping = [
        ["F1.1_E", "F1.2_E"],
        ["F4.3_C", "F4.2_C", "F4.1_C"],
        ]
    
    expected_output = {
        "F1.2_E": ["F1.1_E", "F1.2_E"], 
        "F4.2_C": ["F4.3_C", "F4.2_C", "F4.1_C"],        
    }

    result_json = map_representative_paralog_to_all_redundant_paralogs(
        example_fasta, example_grouping, tmp_path)
    
    assert result_json == expected_output

    example_output_path = Path(tmp_path) / "mapping_redundant_paralogs.json"
    assert example_output_path.exists()

    with open(example_output_path, "r", encoding="utf-8") as f:
        file_contents = json.load(f)

    assert file_contents == expected_output



def test_reduce_seq_collection_to_non_redundancy(tmp_path):
    # mock fasta
    pre_filtered_fasta_content = """>F4.3_C
AAAAAAAAA
>unrelated_header
QQQQQQQQQQ
>F1.1_E
EEEEEEEEEEEE
>unrelated_header2
QQQQQQ
>F1.2_E
EEEEEEEEEEEEEEEEEEE
>F4.1_C
AAAAAAA
>F4.2_C
AAAAAAAAAAAA
>unrelated_header3
QQQQQQQQ
>F5_C
AAA
"""
    pre_filtered_fasta = tmp_path / "pre_filtered.fasta"
    pre_filtered_fasta.write_text(pre_filtered_fasta_content)
    example_grouping = [
        ["F1.1_E", "F1.2_E"],
        ["F4.3_C", "F4.2_C", "F4.1_C"],
        ]
    mapping_json = map_representative_paralog_to_all_redundant_paralogs(
        pre_filtered_fasta, example_grouping, tmp_path)

    post_filtered_fasta_content = """>unrelated_header
QQQQQQQQQQ
>unrelated_header2
QQQQQQ
>F1.2_E
EEEEEEEEEEEEEEEEEEE
>F4.2_C
AAAAAAAAAAAA
>unrelated_header3
QQQQQQQQ
>F5_C
AAA
"""

    expected_filtered_fasta = tmp_path / "post_filtered.fasta"
    expected_filtered_fasta.write_text(post_filtered_fasta_content)

    print(mapping_json)
    reduce_seq_collection_to_non_redundancy(
        input_fasta=pre_filtered_fasta,
        output_dir=tmp_path, 
        redundancy_mapping=mapping_json
    )

    result_output_path = Path(tmp_path) / "paralog_redundancy_filtered.fasta"
    assert result_output_path.exists()

    with open(result_output_path, "r", encoding="utf-8") as f:
        fasta = f.read()

    with open(expected_filtered_fasta, "r", encoding="utf-8") as f:
        expected = f.read()

    assert fasta == expected



# %%
