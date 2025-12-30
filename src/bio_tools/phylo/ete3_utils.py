from ete3 import Tree


def is_fork_node(node:Tree) -> bool:
    """
    Take a ete3 Tree instance (a node) and return true, if node is a forking node


    Background: the ete3 Tree is represented by instances of nodes. 
        Nodes can be categorized into different structural components of the tree:
        - leaf nodes 
        - branch nodes 
            (branch lengths are represented as nodes in ete3. See example below (indicated by **)
        - forking nodes
    e.g.

            /-F7_C
        /-|
    --|    \\-F7.2_B
      |
       *\\-* /-F7_D

    """
    return not node.is_leaf() and len(node.children) >= 2


def fork_distance(leaf: Tree, ancestor: Tree):
    """
    Count number of forking internal nodes between a leaf and an ancestor.
    """
    count = 0
    n = leaf

    while n is not None and n != ancestor:
        n = n.up
        if n is None:
            break
        if is_fork_node(n):
            count += 1

    return count