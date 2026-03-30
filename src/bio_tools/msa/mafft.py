import subprocess
import os
from Bio import AlignIO
from Bio.Seq import Seq

from Bio.Align import MultipleSeqAlignment
from pathlib import Path
import numpy as np

def run_mafft(input_fasta, output_fasta, log_file="mafft.log", keep_log=False):
    """
    Run MAFFT using the command:
        mafft --auto input > output

    Parameters
    ----------
    input_fasta : str
        Path to the input FASTA file.
    output_fasta : str
        Path where the aligned FASTA will be written.
    log_file : str
        Filename to store MAFFT output log.
    keep_log : bool
        If False, the log file will be removed after alignment.

    Returns
    -------
    None
    """
    print("Start mafft")

    # Open log file for writing
    with open(log_file, "w") as lf:
        # Run MAFFT
        result = subprocess.run(
            ["mafft", "--auto", input_fasta],
            stdout=subprocess.PIPE,
            stderr=lf,     # log MAFFT stderr (its progress/messages)
            text=True
        )

    # Write aligned sequences to output file
    with open(output_fasta, "w") as out:
        out.write(result.stdout)

    # Optionally delete log
    if not keep_log and os.path.exists(log_file):
        os.remove(log_file)
    

def trim_msa_by_gap_fraction(input_alignment, output_alignment=None, gap_threshold=0.9):
    """
    Trim an MSA by removing columns with too many gaps.

    Parameters
    ----------
    input_alignment : str or Path
        Path to the input MSA FASTA file.
    output_alignment : str or Path, optional
        Path to save the trimmed alignment. If None, overwrite input.
    gap_threshold : float
        Maximum allowed fraction of gaps per column (0–1).

    Returns
    -------
    trimmed : Bio.Align.MultipleSeqAlignment
        The trimmed alignment.
    """
    print("Start trimming MSA")
    

    input_alignment = Path(input_alignment)
    if output_alignment is None:
        output_alignment = input_alignment
    else:
        output_alignment = Path(output_alignment)

    align = AlignIO.read(input_alignment, "fasta")
    n_seq = len(align)
    aln_len = align.get_alignment_length()

    # Convert to NumPy array (n_seq x aln_len)
    arr = np.array([list(str(rec.seq)) for rec in align], dtype="U1")

    # Compute gap fraction per column
    gap_fraction = np.mean(arr == "-", axis=0)

    # Boolean mask of columns to keep
    keep_mask = gap_fraction <= gap_threshold

    # Slice alignment
    trimmed_arr = arr[:, keep_mask]

    # Build trimmed alignment
    trimmed_records = []
    for rec, row in zip(align, trimmed_arr):
        rec.seq = Seq("".join(row))
        trimmed_records.append(rec)

    trimmed = MultipleSeqAlignment(trimmed_records)
    AlignIO.write(trimmed, output_alignment, "fasta")

    return trimmed