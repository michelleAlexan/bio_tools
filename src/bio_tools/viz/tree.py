from pycirclize import Circos
from ete3 import Tree
from pathlib import Path
from pycirclize import Circos
from bio_tools.phylo.twoODDs import COLORS_2ODD_FUNCTION

from collections import defaultdict
# target gene family is 2ODD but may be any other gene family
GENE_FAM = "2ODD"

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_FOLDER = Path(__file__).resolve().parents[3] / "data"


PATH_BAITS_TREE = DATA_FOLDER / f"{GENE_FAM}s/{GENE_FAM}_char_baits_tree.nwk"

# %%


#-----------find maximal monophyletic subclades ----------
def maximal_monophyletic_clades(tree, target_leaves):
    """
    Return list of clades whose terminals are all in target_leaves,
    but whose parent contains non-target leaves.
    """
    result = []

    for clade in tree.find_clades(order="postorder"):
        terminals = {t.name for t in clade.get_terminals()}
        if not terminals:
            continue

        # clade fully inside target set
        if terminals.issubset(target_leaves):
            parent = tree.root if clade == tree.root else None

            # check if parent is also fully inside (skip if so)
            is_maximal = True
            for parent_candidate in tree.find_clades():
                if clade in parent_candidate.clades:
                    parent = parent_candidate
                    break

            if parent:
                parent_terms = {t.name for t in parent.get_terminals()}
                if parent_terms.issubset(target_leaves):
                    is_maximal = False

            if is_maximal and len(terminals) > 1:
                result.append(terminals)

    return result


# ------------- monophyletic groups with singeltons--------------
def maximal_monophyletic_clades_with_singletons(tree, target_leaves):
    """
    Return list of terminal-name lists.
    Each list is either:
      - a maximal monophyletic clade (>=2 leaves)
      - a singleton leaf not part of any such clade
    """
    clades = []
    covered = set()

    # find maximal monophyletic clades
    for clade in tree.find_clades(order="postorder"):
        terminals = {t.name for t in clade.get_terminals()}
        if not terminals or not terminals.issubset(target_leaves):
            continue

        # check parent
        parent = None
        for parent_candidate in tree.find_clades():
            if clade in parent_candidate.clades:
                parent = parent_candidate
                break

        is_maximal = True
        if parent:
            parent_terms = {t.name for t in parent.get_terminals()}
            if parent_terms.issubset(target_leaves):
                is_maximal = False

        if is_maximal and len(terminals) >= 2:
            clades.append(sorted(terminals))
            covered |= terminals

    # add singletons
    for leaf in target_leaves:
        if leaf not in covered:
            clades.append([leaf])

    return clades

# --------------PLOT and highlight each subclade independently-----------------

def plot_char_2ODD_tree(path_tree: str | Path = PATH_BAITS_TREE):
    """
    Plot highlighted phylogenetic tree. Color to function mapping is done according to the characterized 2ODD tree. 
    """
    circos, tv = Circos.initialize_from_tree(
    path_tree,
    start=10,
    r_lim=(0, 100),
    leaf_label_size=2
    )   

    # --------- prefix -> list of leaf names------------
    function_to_leaves = defaultdict(list)

    for clade in tv.tree.get_terminals():
        if clade.name is None:
            continue
        function = clade.name.split("__")[1]
        function_to_leaves[function].append(clade.name)
    for function, leaves in function_to_leaves.items():
        color = COLORS_2ODD_FUNCTION.get(function)
        if color is None:
            continue

        subclades = maximal_monophyletic_clades_with_singletons(tv.tree, leaves)

        for terminals in subclades:
            tv.highlight(
                list(terminals),
                color=color,
                alpha=0.85,
                lw=1.2,
            )

    fig = circos.plotfig()
    fig.set_dpi(600)

