import subprocess

def run_fasttree(input_alignment: str, output_tree: str, model: str = "-lg"):
    """
    Run FastTree to compute a phylogenetic tree.

    Parameters
    ----------
    input_alignment : str
        Path to the input multiple sequence alignment (FASTA).
    output_tree : str
        Path where the output tree will be written.
    model : str, optional
        Substitution model to use, default is "-lg" (LG model).
    """
    result = subprocess.run(
        ["fasttree", model, input_alignment],
        capture_output=True,
        text=True,
        check=True
    )

    # Take the last line only (the Newick tree)
    tree_line = result.stdout.strip().splitlines()[-1]

    # Write the tree line to the output file
    with open(output_tree, "w") as f:
        f.write(tree_line + "\n")

    return tree_line