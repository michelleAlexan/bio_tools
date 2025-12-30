import pytest
from ete3 import Tree

from bio_tools.phylo.ete3_utils import (
    is_fork_node, 
    fork_distance
    )

def test_is_fork_node():
    tree = Tree("((A:0.3):0.3, (B, C));")
    expected = [
        False, # A
        False, # -- A
        False, # B
        False, # C
               #     /- B
        True,  # --
               #     \\- C

               #         /-----A
               #      --| 
               #        |   /- B
        True   #         \\--
               #            \\- C
    ]
    actual = []
    for n in tree.traverse("postorder"):
        print(n)
        fork_bool = is_fork_node(n)
        actual.append(fork_bool)
    
    assert actual == expected


@pytest.mark.parametrize(
    "input_leaf, expected",  
    [
        (
            "A", 0
        ),
        (
            "B", 1
        )
    ]
)
def test_fork_distance(input_leaf, expected):
    t = Tree("(A:0.3,(B,C));", format=1)
    leaf = t & input_leaf
    result = fork_distance(leaf, t.get_tree_root())
    assert result == expected
