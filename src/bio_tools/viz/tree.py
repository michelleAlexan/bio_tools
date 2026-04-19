#%%
from ete4 import PhyloTree, Tree, NCBITaxa
from ete4.smartview import Layout, TextFace, LegendFace, PropFace
from typing import Callable
from pathlib import Path
from pycirclize import Circos
from collections import defaultdict


def is_char_bait_sequence(leaf_name: str) -> bool:
    return len(leaf_name.split("__")) == 4


VALID_RANKS = {"species", "genus", "family", "order"}
PATH_BAITS_TREE ="/Users/michellealexander/projects/bait_sequence_collection/data/2ODDs/2ODD_char_baits_tree.nwk"

COLORS_CHAR_2ODD_FUNCTION ={
    "GAME31" : "#84cbb6",
    "GAME32" : "#588e82",
    "GAME33" : "#67a790",
    "GAME34" : "#63869c",
    "AOP2" : "#4e7b3a",
    "AOP3" : "#4a7638",
    "DPS" : "#4a7637",
    "GA20ox" : "#6aa84f",
    "C20_GA2ox": "#93c47d",
    "C19_GA2ox": "#b6d7a8",
    "GA2ox" : "#759c63",
    "DAO" : "#9bb78f",
    "GA3ox" : "#b6d7a8",
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
    "E8": "#e06666",
    "GAME40" : "#b05555",
    "D4H" : "#bf7979",
    "BX6" : "#e0bbbb",
    "FNSI" : "#c492cc",
    "FNSI_F3H" : "#c492cc",
    "FNSI_FLS" : "#c492cc",
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
    "FLS_F3H" : "#b4a7d6",
    "DAH": "#784fe1",
    "ANS" : "#8e7cc3",
    "JOX" : "#6fa8dc",
    "ACCO" : "#3d85c6",
    "T6OD" : "#3371a8",
    "COD" : "#316ca2",
    "SRG" : "#316a9f",
    "LBO" : "#2c6190",
    "T2OGD" : "#5D7845",
}


TWO_ODD_COLOR_MAP = {
    "2ODD01": "#c4f4ee",
    "2ODD02": "#34cbc6",
    "2ODD03": "#1b9aa3",
    "2ODD04": "#57c2e0",
    "2ODD05": "#4bdff0",
    "2ODD06": "#90c3c8",
    "2ODD07": "#4f8da8",
    "2ODD08": "#4e9eee",
    "2ODD09": "#0c3fbe",
    "2ODD10": "#D3C8DF",
    "2ODD11": "#8382C4",
    "2ODD11A": "#BAC4F2",
    "2ODD11B": "#10007C",
    "2ODD12": "#555073",
    "2ODD13": "#eaa8e8",
    "2ODD13A": "#c9abc5",
    "2ODD14": "#985C8D",
    "2ODD15": "#da89d1",
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

    "minor_2ODD_cluster": "#999999", 

    "candidate": "#221F1F",
}


_palette = list(TWO_ODD_COLOR_MAP.values())



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

# REGEX = "^(.+?)__(.+?)__(.+?)__(\d+)$" # /r ^(.+?)__(.+?)__(.+?)__(\d+)$
#%%


def _color_for_two_odd_id(two_odd_id: str) -> str:
    """
    Return a colour for the given 2ODD ID.

    - If the ID is in the official colour map, use that colour.
    - Otherwise, assign a deterministic fallback colour drawn from the same palette.
    """
    if two_odd_id in TWO_ODD_COLOR_MAP:
        return TWO_ODD_COLOR_MAP[two_odd_id]

    if not _palette:
        # extremely defensive: if map is somehow empty
        return "#808080"

    # deterministic but “any” colour from the existing palette
    idx = hash(two_odd_id) % len(_palette)
    return _palette[idx]


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
    elif any(x in lineage for x in ["amborellales",
        "nymphaeales",
        "austrobaileyales", 
        "magnoliidae"]):
        return "Basal Angiosperms"
    if "acrogymnospermae" in lineage:
        return "Gymnosperms"
    elif "liliopsida" in lineage:
        return "Monocots"
    elif any(x in lineage for x in ["eudicotyledons", 
                                    "magnoliopsida", 
                                    "mesangiospermae"]):
        return "Dicots"

    else:
        print(f"Plant group couldnt be mapped for node {node.props['sci_name']}")
        print(lineage)


def assign_props_to_leaves(
        tree: Tree, seq_to_two_odd_id: dict | None = None, 
        candidate_headers: set | None = None):

    # --- Assign characterized bait sequence properties ---
    for leaf in tree.leaves():

        if is_char_bait_sequence(leaf.name):
            accession, function, metabolic_pathway, tax_id = leaf.name.split("__")
            leaf.add_props(
                function=function,
                metabolic_pathway=metabolic_pathway
            )

        if seq_to_two_odd_id:
            two_odd_id = seq_to_two_odd_id.get(leaf.name, "None")
            if "minor" in two_odd_id:
                two_odd_id = "minor_2ODD_cluster"
            leaf.add_props(two_odd_id=two_odd_id)
        if candidate_headers and leaf.name in candidate_headers:
            leaf.add_props(two_odd_id="candidate")

        # assign plant group and color based on the plant group
        plant_group = classify_plant(leaf)
        leaf.add_props(
            plant_group=plant_group,
            color=GROUP_COLORS.get(plant_group, None)
        )


def assign_cluster_colors_modern(tree, color_dict=TWO_ODD_COLOR_MAP):
    """
    Assign cluster colors and return an ordered legend mapping.
    """

    # --- First pass: assign colors to leaves ---
    present_ids = set()

    for node in tree.traverse("postorder"):

        if node.is_leaf:
            two_odd_id = node.props.get("two_odd_id")

            if two_odd_id is None:
                color = None
                key = None

            elif "minor" in two_odd_id:
                key = "minor_2ODD_cluster"
                color = color_dict.get(key)

            else:
                key = two_odd_id
                color = color_dict.get(key)

            node.add_props(cluster_color=color)

            if key is not None:
                present_ids.add(key)

            continue

        # --- internal nodes inherit color if uniform ---
        child_colors = {child.props.get("cluster_color") for child in node.children}

        if len(child_colors) == 1:
            color = next(iter(child_colors))
            if color is not None:
                node.add_props(cluster_color=color)

    # --- Second pass: build ordered legend ---
    cluster_colors = {}

    for key in color_dict.keys():  # preserves original order
        if key in present_ids:
            cluster_colors[key] = color_dict[key]

    return cluster_colors
#%%
# t = explorer(PATH_BAITS_TREE, branch_color_mode="function", ultrametric=True, outgroup_leaf="BAW81934__GRS__glucosinolate_biosynthesis__3726")

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
    # Leaf processing
    # --------------------------------------------------
    for leaf in t.leaves():

        # ---- bait detection ----
        if is_char_bait_sequence(leaf.name):
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
            colormap=COLORS_CHAR_2ODD_FUNCTION,
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
            branch_color = COLORS_CHAR_2ODD_FUNCTION.get(func)

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


def explore_2ODD_IDs(tree, show_leaf_dots=True):

    # assign colors per cluster id (cid -> hex color)
    cluster_colors = assign_cluster_colors_modern(tree, TWO_ODD_COLOR_MAP)

    def draw_tree(tree):

        yield LegendFace(
            "Plant groups",
            variable="discrete",
            colormap=GROUP_COLORS
        )

        # only add the cluster legend if we actually have clusters
        if cluster_colors:
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

        # --- remove ALL node dots if disabled ---
        if not show_leaf_dots:
            yield {'dot': {'opacity': 0}}
            
            # still show labels for leaves
            if node.is_leaf:
                yield PropFace('name', position='right')
            return

        # --- internal nodes (when dots enabled) ---
        if not node.is_leaf:
            yield {'dot': {'opacity': 0}}
            return

        # --- leaf dots ---
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
    assign_props_to_leaves(tree)

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
        color = COLORS_CHAR_2ODD_FUNCTION.get(function)
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
