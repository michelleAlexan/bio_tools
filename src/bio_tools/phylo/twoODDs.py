#%%
from bio_tools.utils.io import load_yaml
from pathlib import Path
from Bio import SeqIO
from typing import Iterable

# %%
COLORS_2ODD_FUNCTION ={
    "AOP2" : "#4e7b3a",
    "AOP3" : "#4a7638",
    "DPS" : "#4a7637",
    "GA20ox" : "#6aa84f",
    "C20-GA2ox": "#93c47d",
    "C19-GA2ox": "#b6d7a8",
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
    "LDOX" : "#8e7cc3",
    "JOX" : "#6fa8dc",
    "ACCO" : "#3d85c6",
    "T6OD" : "#3371a8",
    "COD" : "#316ca2",
    "SRG" : "#316a9f",
    "LBO" : "#2c6190",
}


def is_char_bait_sequence(leaf_name: str) -> bool:
    return len(leaf_name.split("__")) == 4

# -------------------- Create iTOL Tree annotation file -------------------- 
# Coloring each sequence of a 2ODD function corresponding to a specified color.
# Addtionally, each characterized 2ODD sequence should be indicated by a star.
def create_2ODD_itol_tree_annotation_files(
    input_source: Path | Iterable[str],
    colorstrip_path: Path = Path("itol_2ODD_colorstrip.txt"),
    symbol_path: Path = Path("itol_2ODD_star_symbols.txt"),
    color_mapping: dict = COLORS_2ODD_FUNCTION,
):
    """
    Create iTOL annotation files for 2ODD phylogenetic trees.

    This function generates:
    1) A DATASET_COLORSTRIP file that colors leaves based on 2ODD function
    2) A DATASET_SYMBOL file that adds a star marker to the same leaves

    Leaf names are expected to follow the format:
        <accession>__<function>__<pathway>__<taxID>

    Parameters
    ----------
    input_source : Path | Iterable[str]
        Either:
        - Path to a FASTA file containing sequences
        - An iterable (list/set) of leaf names

    colorstrip_path : Path
        Output path for the iTOL color strip annotation file

    symbol_path : Path
        Output path for the iTOL symbol (star marker) annotation file

    color_mapping : dict
        Mapping of function name → hex color code

    Notes
    -----
    - Only the first matching function per leaf is used
    - Marker stars are added independently of color strips
    - iTOL does not evaluate logic; all matching is done here
    """

    if not color_mapping:
        raise ValueError("color_mapping must be provided and non-empty")

    # --------------------------------------------------
    # Collect leaf names
    # --------------------------------------------------
    if isinstance(input_source, (list, set, tuple)):
        leaf_names = list(input_source)
    else:
        leaf_names = []
        for rec in SeqIO.parse(str(input_source), "fasta"):
            header = rec.id
            if len(header.split("__")) == 4:
                leaf_names.append(header)

    # --------------------------------------------------
    # Write BRANCH STYLE dataset
    # --------------------------------------------------
    with open(colorstrip_path, "w") as out:
        out.write("DATASET_COLORSTRIP\n")
        out.write("SEPARATOR TAB\n")
        out.write("DATASET_LABEL\t2ODD_Function\n")
        out.write("COLOR\t#000000\n")
        out.write("STRIP_WIDTH\t25\n")
        out.write("MARGIN\t5\n")
        out.write("SHOW_INTERNAL\t0\n")
        out.write("DATA\n")

        for leaf in leaf_names:
            for function, color in color_mapping.items():
                if function in leaf:
                    out.write(f"{leaf}\t{color}\t{function}\n")
                    break


    # --------------------------------------------------
    # Write SYMBOL (STAR) dataset
    # --------------------------------------------------
    with open(symbol_path, "w") as out:
        out.write("DATASET_SYMBOL\n")
        out.write("SEPARATOR COMMA\n")
        out.write("DATASET_LABEL,2ODD_marker\n")
        out.write("COLOR,#000000\n")
        out.write("MAXIMUM_SIZE,5\n")
        out.write("DATA\n")

        for leaf in leaf_names:
            for function in color_mapping:
                if function in leaf:
                    out.write(
                        f"{leaf},3,1,#000000,1,-1\n"
                    )
                    break

path_2ODD_fasta = "/Users/michellealexander/projects/bait_sequence_collection/data/2ODDs/paralog_redundancy_filtered.fasta"
create_2ODD_itol_tree_annotation_files(path_2ODD_fasta)


# %%
