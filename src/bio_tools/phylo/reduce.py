# %%
# python module for reducing (larger) phylogenetic gene trees 
from ete4 import Tree, PhyloTree
from pathlib import Path
from Bio import SeqIO
import json

from bio_tools.files.fasta import filter_fasta

TARGET_RANKS = {"family", "subfamily", "tribe", "genus"}

def reduce_tree_to_family_representatives(tree: Tree | PhyloTree, fasta_path: str, output_dir: str):
    """
    Collapse clades ONLY at ranks:
    family, subfamily, tribe, genus

    For each such clade:
    - collect all descendant leaves
    - find the longest sequence
    - map representative -> all members

    save mapping and filtered fasta file to output dir as json and fasta, respectively.

    Return mapping, output_fasta_path, output_json_path
    """

    records = list(SeqIO.parse(fasta_path, "fasta"))
    record_dict = {r.id: str(r.seq) for r in records}

    clades_to_collapse = []
    visited = set()

    for node in tree.traverse("levelorder"):
        rank = node.props.get("rank")

        if rank == "clade" or rank is None:
            continue

        if rank in TARGET_RANKS:

            # skip if already covered by a higher-level collapse
            if node in visited:
                continue

            leaves = list(node.leaf_names())
            clades_to_collapse.append(leaves)

            # mark entire subtree as visited
            for desc in node.traverse():
                visited.add(desc)

    # Build mapping: longest_sequence -> all_clade_members
    rep_to_members = {}

    for clade in clades_to_collapse:
        longest_record = None
        max_len = -1

        for leaf in clade:
            aa_seq = record_dict.get(leaf)
            if aa_seq is None:
                continue  # avoid crash if missing

            if len(aa_seq) > max_len:
                max_len = len(aa_seq)
                longest_record = leaf

        if longest_record is not None:
            rep_to_members[longest_record] = clade

    # save mapping and filtered fasta file to output dir as json and fasta, respectively.
    output_fasta_path = Path(output_dir) / "reduced_fam_level.fasta"
    output_json_path = Path(output_dir) / "reduced_fam_level.json"

    # seqs to discard

    seqs_to_discard = []
    for rep, members in rep_to_members.items():
        for member in members:
            if member != rep:
                seqs_to_discard.append(member)
    
    filter_fasta(
        input_fasta=fasta_path,
        output_fasta=output_fasta_path,
        accessions=seqs_to_discard,
        mode="remove"
    )
        
    with open(output_json_path, "w") as f:
        json.dump(rep_to_members, f, indent=4)
    

    return rep_to_members, output_fasta_path, output_json_path


