#%%
from pycirclize import Circos
from ete4 import PhyloTree
from ete4.smartview import Layout, TextFace, LegendFace
from typing import Callable
from pathlib import Path
from pycirclize import Circos
from bio_tools.phylo.twoODDs import COLORS_2ODD_FUNCTION

from ete4.smartview import Layout, TextFace, PropFace, BASIC_LAYOUT

from collections import defaultdict

PATH_BAITS_TREE ="/Users/michellealexander/projects/bait_sequence_collection/data/2ODDs/2ODD_char_baits_tree.nwk"


GROUP_COLORS = {
    "Algae": "#574104",
    "Lycophytes": "#ab730c",
    "Liverworts": "#faaf00",
    "Ferns": "#c1d717",
    "Mosses": "#89be86",
    "Gymnosperms": "#076247",
    "Early Angiosperms": "#3BA0BC",
    "Monocots": "#7502d9",
    "Dicots": "#edc5ec",
}

COLORS_2ODD_FUNCTION ={
    "AOP2" : "#438195",
    "AOP3" : "#438195",
    "DPS" : "#377672",
    "GA20ox" : "#6aa84f",
    "C20_GA2ox": "#406732",
    "C19_GA2ox": "#b6d7a8",
    "DAO" : "#539185",
    "GA3ox" : "#64b541",
    "GA13ox" : "#a0bd94",
    "GA7ox" : "#e8eed0",
    "2ODDC23" : "#fff2cc",
    "LFS" : "#f4e8c3",
    "2OG1" : "#ffe599",
    "C2H" : "#ffd966",
    "F6H" : "#ec8612",
    "S8H" : "#f57829",
    "GSLOH" : "#ec640f",
    "GRS" : "#e5ac00",
    "TIIAS" : "#af8300",
    "D4H" : "#bf7979",
    "BX6" : "#e0bbbb",
    "FNSI" : "#c492cc",
    "F3H" : "#f4cccc",
    "H6H" : "#e9d0db",
    "IDS" : "#cd87a6",
    "SLC" : "#cfe2f3",
    "GIM" : "#d0d2e5",
    "M2H" : "#c27ba0",
    "M2H_weak" : "#e688b8",
    "DMR6" : "#c2d3e2",
    "S5H" : "#abbbc9",
    "S3H" : "#95a3af",
    "FLS" : "#b4a7d6",
    "DAH": "#784fe1",
    "ANS" : "#8e7cc3",
    "JOX" : "#6fa8dc",
    "ACCO" : "#3d85c6",
    "T6OD" : "#3371a8",
    "COD" : "#316ca2",
    "SRG" : "#316a9f",
    "LBO" : "#2c6190",
    "NCS" : "#1a3852",

}

# REGEX = "^(.+?)__(.+?)__(.+?)__(\d+)$" # /r ^(.+?)__(.+?)__(.+?)__(\d+)$
def is_char_bait_sequence(leaf_name: str) -> bool:
    return len(leaf_name.split("__")) == 4


def classify_plant(node):
    lineage = [t.lower() for t in node.props["named_lineage"]]

    if "zygnematophyceae" in lineage:
        return "Algae"
    elif "lycopodiopsida" in lineage:
        return "Lycophytes"    #Non-seed vascular plants
    elif "polypodiopsida" in lineage:
        return "Ferns"
    elif "bryophyta" in lineage or "anthocerotophyta" in lineage:
        return "Mosses"
    elif "marchantiophyta" in lineage:
        return "Liverworts"
    if "acrogymnospermae" in lineage:
        return "Gymnosperms"
    elif "liliopsida" in lineage:
        return "Monocots"
    elif any(x in lineage for x in ["eudicotyledons", 
                                    "magnoliopsida", 
                                    "mesangiospermae"]):
        return "Dicots"


    elif any(x in lineage for x in ["amborellales",
        "nymphaeales",
        "austrobaileyales", 
        "magnoliidae"]):
        return "Basal Angiosperms"
    else:
        print(f"Plant group couldnt be mapped for node {node.props["sci_name"]}")
        print(lineage)


#%%



def explorer(
    newick: str|Path,
    branch_color_mode: str = "plant",  # "plant" or "function"
    sp_naming_function: Callable[[str], str] = lambda name: name.split("__")[-1],
    ladderize: bool = True,
    ultrametric: bool = False, 
    outgroup_leaf: str = None
):
    """
    Explorer for phylogenetic trees with flexible clade coloring.

    Parameters
    ----------
    newick : str | Path
        Path to Newick file or Newick string.
    branch_color_mode : str
        "plant"     → branches colored by plant class
        "function"  → branches colored by 2ODD function
    sp_naming_function : callable
        Function extracting taxid from leaf name.
    ladderize : bool
        Ladderize tree before visualization.
    """

    if branch_color_mode not in {"plant", "function"}:
        raise ValueError("branch_color_mode must be 'plant' or 'function'")

    # --------------------------------------------------
    # Load tree
    # --------------------------------------------------
    if isinstance(newick, (str, Path)) and Path(str(newick)).exists():
        t = PhyloTree(open(newick), sp_naming_function=sp_naming_function)
    elif isinstance(newick, str):
        t = PhyloTree(newick, sp_naming_function=sp_naming_function)
    else:
        raise ValueError("newick must be valid path or Newick string")

    # --------------------------------------------------
    # Taxonomy annotation
    # --------------------------------------------------
    try:
        t.annotate_ncbi_taxa(taxid_attr="species")
    except Exception as e:
        print(f"[WARNING] Taxonomy annotation failed: {e}")

    if ladderize:
        t.ladderize()

    if outgroup_leaf:
        t.set_outgroup(node=t[outgroup_leaf])
    
    if ultrametric:
        t.to_ultrametric()

    # --------------------------------------------------
    # Detect characterized bait sequence
    # --------------------------------------------------
    def is_char_bait(name: str):
        return len(name.split("__")) == 4

    # --------------------------------------------------
    # Leaf processing
    # --------------------------------------------------
    for leaf in t.leaves():

        # ---- bait detection ----
        if is_char_bait(leaf.name):
            try:
                accession, function, pathway, taxid = leaf.name.split("__")
                leaf.add_props(
                    function=function,
                    metabolic_pathway=pathway,
                    is_bait=True,
                )
            except ValueError:
                leaf.add_props(is_bait=False)
        else:
            leaf.add_props(is_bait=False)

        # ---- plant classification ----
        try:
            plant_group = classify_plant(leaf)
        except Exception:
            plant_group = None

        leaf.add_props(
            plant_group=plant_group,
            plant_color=GROUP_COLORS.get(plant_group),
        )

    # --------------------------------------------------
    # Propagate function annotation to internal nodes
    # (only if leaves contain function annotation)
    # --------------------------------------------------
    if branch_color_mode == "function":
        for node in t.traverse("postorder"):
            if node.is_leaf:
                continue

            child_functions = {
                child.props.get("function")
                for child in node.leaves()
                if child.props.get("function") is not None
            }

            if len(child_functions) == 1:
                node.add_props(function=child_functions.pop())

    # --------------------------------------------------
    # Layout
    # --------------------------------------------------
    def draw_tree(tree):

        # ---- Plant legend ----
        yield LegendFace(
            "Plant groups",
            variable="discrete",
            colormap=GROUP_COLORS,
        )

        # ---- Function legend ----
        yield LegendFace(
            "2ODD Functions",
            variable="discrete",
            colormap=COLORS_2ODD_FUNCTION,
        )

        yield {"node-height-min": 1.0, 
               "shape": "circular",}

    def draw_node(node):

        # --------------------------------------------------
        # Branch coloring
        # --------------------------------------------------
        if branch_color_mode == "plant":
            branch_color = node.props.get("plant_color")

        elif branch_color_mode == "function":
            func = node.props.get("function")
            branch_color = COLORS_2ODD_FUNCTION.get(func)

        else:
            branch_color = None

        if branch_color:
            yield {
                "hz-line": {
                    "stroke": branch_color,
                    "stroke-width": 3,
                },
                "vt-line": {
                    "stroke": branch_color,
                    "stroke-width": 3,
                },
            }

        # --------------------------------------------------
        # Remove internal node dots
        # --------------------------------------------------
        if not node.is_leaf:
            yield {"dot": {"opacity": 0}}
            return

        # --------------------------------------------------
        # Leaf styling
        # Yellow square with plant color fill
        # --------------------------------------------------
        plant_color = node.props.get("plant_color")

        yield {
            "dot": {
                "shape": "square",
                "radius": 7,
                "fill": plant_color or "white",
                # "stroke": plant_color,
                # "stroke-width": 3,
                "opacity": 1,
                "fill-opacity": 1,

            }
        }

        yield PropFace("name", position="right")

    layout = Layout(
        "Flexible clade coloring",
        draw_tree=draw_tree,
        draw_node=draw_node,
    )

    # --------------------------------------------------
    # Launch explorer
    # --------------------------------------------------
    t.explore(
        layouts=[layout],
        show_leaf_name=False,
        show_popup_props=[
            "name",
            "sci_name",
            "taxid",
            "rank",
            "plant_group",
            "function",
            "metabolic_pathway",
        ],
    )

    return t

#%%
t = explorer(PATH_BAITS_TREE, branch_color_mode="function", ultrametric=True, outgroup_leaf="BAW81934__GRS__glucosinolate_biosynthesis__3726")

#%%
accessions = list(t[1,1,1,1,1,1,0,1,1,1,1,0].leaves())































#%%

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



#%%
