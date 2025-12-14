import subprocess
from pathlib import Path

def run_fasttree(
    input_alignment: Path,
    output_tree: Path,
    model: str = "-lg",
    fasttree_bin: str = "fasttree",
):
    """
    Run fasttree to compute a phylogenetic tree.

    Parameters
    ----------
    input_alignment : Path
        Path to aligned FASTA file.
    output_tree : Path
        Path to write Newick tree.
    model : str
        Substitution model (e.g. "-lg", "-gtr").
    fasttree_bin : str
        fasttree executable name.
    """

    input_alignment = Path(input_alignment)
    output_tree = Path(output_tree)

    if not input_alignment.exists():
        raise FileNotFoundError(input_alignment)

    cmd = [
        fasttree_bin,
        model,
        str(input_alignment),
    ]

    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,   # tree is written here
        stderr=subprocess.PIPE,   # diagnostics here
        text=True,
    )

    if result.returncode != 0:
        raise RuntimeError(
            f"fasttree failed with exit code {result.returncode}\n"
            f"STDERR:\n{result.stderr}"
        )

    # FastTree prints exactly one Newick tree
    tree = result.stdout.strip()

    if not tree.startswith("("):
        raise ValueError(
            "fasttree output does not look like a Newick tree:\n"
            f"{tree[:200]}"
        )

    output_tree.write_text(tree + "\n")

