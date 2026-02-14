#%%
from pycirclize import Circos
from ete4 import Tree, PhyloTree
from ete4.smartview import Layout, TextFace, LegendFace

from pathlib import Path
from pycirclize import Circos
from bio_tools.phylo.twoODDs import COLORS_2ODD_FUNCTION
import re
from math import pi
from ete4.smartview import Layout, TextFace, BASIC_LAYOUT

from collections import defaultdict
# target gene family is 2ODD but may be any other gene family
GENE_FAM = "2ODD"

PACKAGE_ROOT = Path(__file__).resolve().parents[1]
DATA_FOLDER = Path(__file__).resolve().parents[3] / "data"


PATH_BAITS_TREE = DATA_FOLDER / f"{GENE_FAM}s/{GENE_FAM}_char_baits_tree.nwk"


GROUP_COLORS = {
    "Lycophytes": "#f6a50d",
    "Liverworts": "#a6761d",
    "Ferns": "#c1d717",
    "Mosses": "#89be86",
    "Gymnosperms": "#076247",
    "Early Angiosperms": "#3BA0BC",
    "Monocots": "#7502d9",
    "Dicots": "#edc5ec",
    "Other": "#5F5F5F"
}

def is_char_bait_sequence(leaf_name: str) -> bool:
    return len(leaf_name.split("__")) == 4


def classify_plant(node):
    lineage = [t.lower() for t in node.props["named_lineage"]]

    if "lycopodiopsida" in lineage:
        return "Lycophytes"    #Non-seed vascular plants
    elif "polypodiopsida" in lineage:
        return "Ferns"
    elif "bryophyta" in lineage:
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
        return "Other"

#%%

TREE_PATH = Path("/Users/michellealexander/Documents/ingroup_outgroup/B3/post_B3_tree.nwk")
t = PhyloTree(open(TREE_PATH), sp_naming_function=lambda name: name.split('__')[-1])
tax2names, tax2lineages, tax2rank = t.annotate_ncbi_taxa(taxid_attr='species')

#%%
for leaf in t.leaves():
    if is_char_bait_sequence(leaf.name):
        accession, function, metabolic_pathway, tax_id = leaf.name.split("__")
        leaf.add_props(function= function, 
                       metabolic_pathway= metabolic_pathway)
        

    plant_group = classify_plant(leaf)
    leaf.add_props(plant_group=plant_group, 
                    color=GROUP_COLORS.get(plant_group, "none"))
    print(leaf.props["color"])


def draw_tree(tree):
    yield LegendFace(
        "Plant groups",
        variable="discrete",
        colormap=GROUP_COLORS
    )


def draw_node(node):

    color = node.props.get("color")

    if color:
        yield {'box': {'fill': color, 'opacity': 0.6}}

    if node.is_leaf:
        yield {
            'dot': {
                'shape': 'circle',
                'radius': 8,
                'fill': color,
                'stroke': color
            }
        }


layout = Layout(
    "Plant group layout",
    draw_tree=draw_tree,
    draw_node=draw_node
)

t.explore(
    layouts=[layout],
    show_popup_props=[
        'name', 'sci_name', 'taxid',
        'rank', 'plant_group'
    ]
)

#%%

#print(t.to_str(props=['name', 'sci_name', 'taxid', 'rank']))



#%%


















































#%%

pattern = re.compile(r'^(.+?)__(.+?)__(.+?)__(\d+)$')

# ---- TREE STYLE (GLOBAL) ----
tree_style = {
    'shape': 'circular',
    'radius': 0,
    'angle-start':  -5 * pi / 6     # where the tree starts
}

# ---- NODE STYLING ----
def draw_node(node):
    if node.is_leaf:

        # Always draw the name
        yield TextFace(node.name,
                       fs_min=10, fs_max=14,  # font size
                       style={"fill": "black"},
                       position="aligned")  # aligns nicely next to the node

        # Optional decorations
        if pattern.match(node.name):
            yield {"box": {"fill": "red"}}
            yield TextFace("★", fs_min=12, fs_max=16, style={"fill": "orange"}, position="aligned")

layout = Layout(
    name="highlight_pattern",
    draw_tree=tree_style,
    draw_node=draw_node
)

# load your tree
TREE_PATH = Path("/Users/michellealexander/Documents/ingroup_outgroup/B3/post_B3_tree.nwk")
t = Tree(open(TREE_PATH))

t.explore(layouts=[layout])
#%%

#%%
def taxon_block_layout(node):
    dist = node.get_distance(t)
    if dist < 0.3:
        rank_to_show = "genus"
    else:
        rank_to_show = "order"

    name = getattr(node, rank_to_show, None)
    if name:
        color = order_colors.get(name, "#CCCCCC")
        node.add_face(RectFace(15,10, color, color), column=1, position="branch-right")
        if node.is_leaf():
            node.add_face(TextFace(node.sci_name, fsize=8), column=2, position="aligned")

from ete4 import BarChartFace

# precompute counts per node (example: counts per order)
for node in t.traverse():
    counts = Counter(ch.order for ch in node.iter_leaves())
    node.add_feature("order_counts", counts)

def heatmap_layout(node):
    if hasattr(node, "order_counts"):
        labels = list(node.order_counts.keys())
        values = list(node.order_counts.values())
        face = BarChartFace(values, labels=labels, width=80, height=10, colors=[order_colors[l] for l in labels])
        node.add_face(face, column=3, position="branch-right")








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
