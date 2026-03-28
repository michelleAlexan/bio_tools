#%%
from ete4 import PhyloTree, Tree, NCBITaxa
from ete4.smartview import Layout, TextFace, LegendFace, PropFace
from typing import Callable
from pathlib import Path
from pycirclize import Circos
from bio_tools.phylo.twoODDs import COLORS_2ODD_FUNCTION
from collections import defaultdict


VALID_RANKS = {"species", "genus", "family", "order"}
PATH_BAITS_TREE ="/Users/michellealexander/projects/bait_sequence_collection/data/2ODDs/2ODD_char_baits_tree.nwk"



COL_2ODD_CLADES = {
    "2ODD01": "#c4f4ee",
    "2ODD02": "#34cbc6",
    "2ODD03": "#1b9aa3",
    "2ODD04": "#57c2e0",
    "2ODD05": "#4bdff0",
    "2ODD06": "#90c3c8",
    "2ODD07": "#4f8da8",
    "2ODD08": "#4e9eee",
    "2ODD09": "#0c3fbe",
    "2ODD10": "#10007C",

    "2ODD11": "#BAC4F2",
    "2ODD11A": "#8382C4",
    "2ODD11B": "#555073",

    "2ODD12": "#D3C8DF",
    "2ODD13": "#eaa8e8",

    "2ODD14": "#c9abc5",
    "2ODD14A": "#da89d1",

    "2ODD15": "#985C8D",
    "2ODD16": "#682c69",
    "2ODD17": "#985fc9",
    "2ODD18": "#5d3c79",
    "2ODD19": "#db1dc2",
    "2ODD20": "#ab2599",
    "2ODD21": "#df3960",
    "2ODD22": "#794058",
    "2ODD23": "#48353D",
    "2ODD24": "#A91010",
    "2ODD25": "#4F0303",
    "2ODD26": "#553CC6",
    "2ODD27": "#3CC6AF",
    "2ODD28": "#B4E1D6",
    "2ODD29": "#71D3AF",
    "2ODD30": "#137549",
    "2ODD31": "#143E1A",
    "2ODD32": "#6b9113",
    "2ODD33": "#218e05",
    "2ODD34": "#26ba09",
    "2ODD35": "#cce066",
    "2ODD36": "#c1c80c",
    "2ODD37": "#ffd966",
    "2ODD38": "#dcbf6f",
    "2ODD39": "#a4830b",
    "2ODD40": "#d67906",
    "2ODD41": "#e33c09",

    "2ODD_minor_clades": "#999999", 

    "candidate": "#010101",
}


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
        print(f"Plant group couldnt be mapped for node {node.props['sci_name']}")
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



def assign_props_to_leaves(tree: Tree, seq_to_two_odd_id: dict | None = None):

    # --- Assign characterized bait sequence properties ---
    for leaf in tree.leaves():

        if is_char_bait_sequence(leaf.name):
            accession, function, metabolic_pathway, tax_id = leaf.name.split("__")
            leaf.add_props(
                function=function,
                metabolic_pathway=metabolic_pathway
            )

        if seq_to_two_odd_id:
            two_odd_id = seq_to_two_odd_id.get(leaf.name, "minor_2ODD_cluster")
            leaf.add_props(two_odd_id=two_odd_id)

        # assign plant group and color based on the plant group
        plant_group = classify_plant(leaf)
        leaf.add_props(
            plant_group=plant_group,
            color=GROUP_COLORS.get(plant_group, None)
        )


def explore_tree_plant_groups(tree: Tree):

    # --- Layout functions ---
    def draw_tree(tree):
        yield LegendFace(
            "Plant groups",
            variable="discrete",
            colormap=GROUP_COLORS
        )
        
        yield {'node-height-min': 1.0}


    def draw_node(node):

        color = node.props.get("color")

        # --- Background color for branches ---
        if color:
            yield {
                'box': {
                    'fill': color,
                    'opacity': 0.6
                }
            }

        # --- Remove internal node dots ---
        if not node.is_leaf:
            yield {'dot': {'opacity': 0}}
            return

        # --- Leaf styling ---
        yield {
            'dot': {
                'shape': 'circle',
                'radius': 6,
                'fill': color,
                'stroke': color
            }
        }
        yield PropFace('name', position='right')



    layout = Layout(
        "Plant group layout",
        draw_tree=draw_tree,
        draw_node=draw_node
    )


    tree.explore(
        layouts=[layout],
        show_leaf_name=True,   
        show_popup_props=[
            'name', 'sci_name', 'two_odd_id', 'taxid', 'named_lineage',
            'rank', 'plant_group'
        ]
    )




def load_treecluster_assignments(tree: Tree, cluster_file: str):

    cluster_map = {}

    with open(cluster_file) as f:
        for line in f:

            if line.startswith("SequenceName"):
                continue

            name, cluster = line.strip().split()
            cluster_map[name] = int(cluster)

    for leaf in tree.leaves():
        two_odd_id = cluster_map.get(leaf.name, -1)
        leaf.add_props(two_odd_id=two_odd_id)

def assign_cluster_colors_modern(tree, color_dict):
    """
    Assign cluster colors using predefined COL_2ODD_CLADES dictionary.
    """
    clusters = {}
    for leaf in tree.leaves():
        cid = leaf.props.get("two_odd_id", -1)
        if cid == -1:
            continue
        clusters.setdefault(cid, []).append(leaf)

    cluster_colors = {}

    # Prepare ordered color keys (exclude minor clades for now)
    color_keys = [k for k in color_dict.keys() if k != "2ODD_minor_clades"]
    color_keys_sorted = sorted(color_keys)

    for i, cid in enumerate(sorted(clusters)):
        if i < len(color_keys_sorted):
            key = color_keys_sorted[i]
        else:
            # fallback if more clusters than colors
            key = "2ODD_minor_clades"

        cluster_colors[cid] = color_dict[key]

    # propagate cluster color upward only when all descendants belong to same cluster
    for node in tree.traverse("postorder"):
        if node.is_leaf:
            cid = node.props.get("two_odd_id", -1)
            if cid != -1:
                node.add_props(cluster_color=cluster_colors[cid])
            continue

        child_colors = {child.props.get("cluster_color") for child in node.children}

        if len(child_colors) == 1:
            color = next(iter(child_colors))
            if color:
                node.add_props(cluster_color=color)

    return cluster_colors



def explore_tree_cluster_clades(tree):

    cluster_colors = assign_cluster_colors_modern(tree, COL_2ODD_CLADES)

    def draw_tree(tree):

        yield LegendFace(
            "Plant groups",
            variable="discrete",
            colormap=GROUP_COLORS
        )

        yield LegendFace(
            "TreeCluster clades",
            variable="discrete",
            colormap=cluster_colors
        )

        yield {'node-height-min': 1.0}

    def draw_node(node):

        cluster_color = node.props.get("cluster_color")

        # --- color branches by cluster ---
        if cluster_color:
            yield {
                'hz-line': {'stroke': cluster_color, 'stroke-width': 3},
                'vt-line': {'stroke': cluster_color, 'stroke-width': 3}
            }

        # --- internal nodes ---
        if not node.is_leaf:
            yield {'dot': {'opacity': 0}}
            return

        # --- leaf dots = plant group ---
        plant_color = node.props.get("color")

        yield {
            'dot': {
                'shape': 'circle',
                'radius': 6,
                'fill': plant_color,
                'stroke': plant_color
            }
        }

        yield PropFace('name', position='right')

    layout = Layout(
        "Cluster clade layout",
        draw_tree=draw_tree,
        draw_node=draw_node
    )

    tree.explore(
        layouts=[layout],
        show_popup_props=[
            'name',
            'two_odd_id',
            'plant_group',
            'sci_name',
            'taxid'
        ]
    )



# %%
def assign_plant_group_props(tree:Tree |PhyloTree):
    """
    Take a ete4 Tree / PhyloTree and assign plant group properties to the leaves based on their taxonomic lineage.
        - Algae
        - Lycophytes (non-seed vascular plants)
        - Liverworts
        - Mosses
        - Ferns
        - Gymnosperms
        - Basal Angiosperms
        - Monocots
        - Dicots (eudicots)
    Note that the tree.annotate_ncbi_taxa() function must have been run beforehand to populate the "named_lineage" property for each leaf, 
    which contains the taxonomic lineage as a list of taxonomic names. 
    """

    GROUP_COLORS = {
            "Algae": "#574104",
            "Lycophytes": "#ab730c",
            "Liverworts": "#faaf00",
            "Ferns": "#c1d717",
            "Mosses": "#89be86",
            "Gymnosperms": "#076247",
            "Basal Angiosperms": "#3BA0BC",
            "Monocots": "#7502d9",
            "Dicots": "#edc5ec",
        }
    
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
        elif any(x in lineage for x in ["amborellales",
            "nymphaeales",
            "austrobaileyales", 
            "magnoliidae"]):
            return "Basal Angiosperms"
        elif "liliopsida" in lineage:
            return "Monocots"
        elif any(x in lineage for x in ["eudicotyledons", 
                                        "magnoliopsida", 
                                        "mesangiospermae"]):
            return "Dicots"


        else:
            print(f"Plant group couldnt be mapped for node {node.props['sci_name']}")
            print(lineage)

    for leaf in tree.leaves():
        plant_group = classify_plant(leaf)
        leaf.add_props(
            plant_group=plant_group,
            color=GROUP_COLORS.get(plant_group, None)
        )



def build_rank_level_tree(
    taxids: set[int],
    rank: str = "species",
    clean_names: bool = True
) -> Tree:
    """
    Build a taxonomy tree at a specified rank and annotate leaves with
    the number of species represented.

    Example leaf name:
        Poales (13)

    Parameters
    ----------
    taxids : set[int]
        Input NCBI taxonomic IDs

    rank : str
        One of {"species", "genus", "family", "order"}

    clean_names : bool
        Whether to sanitize node names

    Returns
    -------
    Tree
        Rank-level taxonomy tree with counts
    """

    if rank not in VALID_RANKS:
        raise ValueError(f"rank must be one of {VALID_RANKS}")

    ncbi = NCBITaxa()

    # --------------------------------------------------
    # Step 1: lineage + rank info
    # --------------------------------------------------
    lineage_dict = ncbi.get_lineage_translator(list(taxids))

    all_lineage_taxids = {
        tid for lineage in lineage_dict.values() for tid in lineage
    }

    rank_dict = ncbi.get_rank(all_lineage_taxids)
    name_dict = ncbi.get_taxid_translator(all_lineage_taxids)

    # --------------------------------------------------
    # Step 2: map each input → target rank taxid
    # AND count how many species map to each
    # --------------------------------------------------
    projected_taxids = set()
    taxid_counts = defaultdict(int)

    for taxid, lineage in lineage_dict.items():

        target_taxid = None

        for t in lineage:
            if rank_dict.get(t) == rank:
                target_taxid = t
                break

        if target_taxid is None:
            target_taxid = taxid  # fallback

        projected_taxids.add(target_taxid)
        taxid_counts[target_taxid] += 1  # count species

    # --------------------------------------------------
    # Step 3: build topology
    # --------------------------------------------------
    tree = ncbi.get_topology(list(projected_taxids))

    # attach taxids
    for leaf in tree.leaves():
        leaf_taxid = int(leaf.name)
        leaf.add_props(taxid=leaf_taxid)

    # annotate taxonomy
    tree.annotate_ncbi_taxa(taxid_attr="taxid")

    # --------------------------------------------------
    # Step 4: rename leaves with counts
    # --------------------------------------------------
    for leaf in tree.leaves():
        taxid = leaf.props["taxid"]
        sci_name = leaf.props.get("sci_name", str(taxid))

        count = taxid_counts.get(taxid, 1)

        if clean_names:
            sci_name = (
                sci_name
                .replace(":", "_")
                .replace("(", "")
                .replace(")", "")
                .replace(" ", "_")
            )

        if rank == "species":
            leaf.name = sci_name
        else:
            leaf.name = f"{sci_name} ({count})"

    # --------------------------------------------------
    # Step 5: plant group annotation
    # --------------------------------------------------
    assign_plant_group_props(tree)

    # --------------------------------------------------
    # Step 6: formatting
    # --------------------------------------------------
    tree.ladderize()
    tree.to_ultrametric()

    return tree


def explore_taxonomy_tree(tree):

    def draw_node(node):

        if not node.is_leaf:
            return

        color = node.props.get("color", "#cccccc")

        yield {
            'dot': {
                'shape': 'circle',
                'radius': 6,
                'fill': color,
                'stroke': color
            }
        }

        yield PropFace('name', position='right')

    layout = Layout("Plant groups", draw_node=draw_node)

    tree.explore(layouts=[layout])


# Define your plant group colors
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




def make_smartview_layout(tree):
    """
    Layout for ete4 smartview render_sm:
    - Internal nodes invisible (dots), branches remain
    - Leaf dots colored by plant group
    - Leaf names displayed
    - Legend for plant groups
    """

    # --- Tree-wide style / legend ---
    def draw_tree(tree):
        legend = LegendFace(
            "Plant groups",
            variable="discrete",
            colormap=GROUP_COLORS  # your dict mapping group -> color
        )
        # Optional: change font size
        legend.label_size = 12
        legend.title_size = 14

        # manually position
        legend.x = 1.0
        legend.y = 0.0
        legend.anchor = "top-right"

        # Branch style (applied globally)
        tree_style = {
            'hz-line': {'stroke': '#000', 'stroke-width': 1},
            'vt-line': {'stroke': '#000', 'stroke-width': 1},
            'node-height-min': 1.0
        }

        return [legend, tree_style]

    # --- Node-specific style ---
    def draw_node(node):
        if node.is_leaf:
            style = {
                'dot': {
                    'shape': 'circle',
                    'radius': 6,
                    'fill': node.props.get('color', '#000'),
                    'stroke': node.props.get('color', '#000')
                }
            }
            # Return dot + leaf name
            return [style, PropFace('name', position='right')]
        else:
            # Internal node dots invisible, branches remain
            return {'dot': {'opacity': 0, 'radius': 0}}

    return Layout(
        name="Plant group layout",
        draw_tree=draw_tree,
        draw_node=draw_node
    )



def render_tree_sm(tree, output_path):
    layout = make_smartview_layout(tree)

    n_leaves = len(list(tree.leaves()))
    height_px = max(2000, n_leaves * 1)

    tree.render_sm(
        output_path,
        layouts=[layout],
        w=1500,
        h=height_px)
#%%
# t = explorer(PATH_BAITS_TREE, branch_color_mode="function", ultrametric=True, outgroup_leaf="BAW81934__GRS__glucosinolate_biosynthesis__3726")

#%%
#accessions = list(t[1,1,1,1,1,1,0,1,1,1,1,0].leaves())


#%% REDUCED PARALOGOUS SEQUENCES

path_red_tree = "/Users/michellealexander/projects/bait_sequence_collection/data/2ODDs/FILTERED_2ODD_char_baits_tree.nwk"
# t = explorer(path_red_tree, branch_color_mode="function", ultrametric=True, outgroup_leaf="BAW81934__GRS__glucosinolate_biosynthesis__3726")




#%%




























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
