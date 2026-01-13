import pandas as pd
from Bio import SeqIO
from typing import Literal
from pathlib import Path
import json
import re

def filter_fasta(input_fasta, output_fasta, accessions, mode="remove"):
    """
    Filter sequences in a FASTA file based on a set of accession strings.

    Parameters
    ----------
    input_fasta : str
        Path to input FASTA.
    output_fasta : str
        Path to output FASTA.
    accessions : set
        Set of accession strings to match in headers.
    mode : str
        "remove" = drop sequences containing any accession (default)
        "keep"   = keep ONLY sequences containing any accession

    """
    mode = mode.lower()
    if mode not in {"remove", "keep"}:
        raise ValueError("mode must be 'remove' or 'keep'")

    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            header = record.description
            match = any(acc in header for acc in accessions)

            if mode == "remove" and match:
                continue
            if mode == "keep" and not match:
                continue

            SeqIO.write(record, out, "fasta")

def find_duplicate_headers(fasta_path):
    seen = set()
    duplicates = []

    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description
        if header in seen:
            duplicates.append(header)
        else:
            seen.add(header)

    return duplicates


def df_to_fasta(
        df: pd.DataFrame, 
        path: str, 
        header_cols=None,
        seq_col: str = "aa_sequence",
        wrap: int = 80
    ):
    """
    Export a DataFrame to FASTA format.

    Parameters
    ----------
    df : pd.DataFrame
        Input table containing sequence and metadata.
    path : str
        Output FASTA file.
    header_cols : list of str
        Columns to include in the FASTA header in order.
    seq_col : str
        Column containing the amino acid sequence.
    wrap : int
        Number of characters per FASTA line.
    """
    if header_cols is None:
        header_cols = ["accession", "function", "metabolic_function", "organism"]

    with open(path, "w") as f:
        for _, row in df.iterrows():

            # Build header from selected columns
            header_parts = []
            for col in header_cols:
                val = str(row.get(col, "NA")).replace(" ", "_")
                header_parts.append(val)

            header = ">" + "__".join(header_parts)

            # Extract sequence
            seq = str(row[seq_col]).strip()

            # wrap sequence lines
            wrapped = "\n".join(seq[i:i+wrap] for i in range(0, len(seq), wrap))

            f.write(f"{header}\n{wrapped}\n")


def fasta_to_df(path: str) -> pd.DataFrame:
    """
    Reads a FASTA file where headers are:
    >accession$function$metabolic_function$organism
    and returns a DataFrame.
    """
    records = []
    accession, function, metabolic_function, organism = None, None, None, None
    seq_lines = []

    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                # save previous record
                if accession is not None:
                    records.append({
                        "accession": accession,
                        "function": function,
                        "metabolic_function": metabolic_function,
                        "organism": organism,
                        "aa_sequence": "".join(seq_lines)
                    })
                # parse new header
                header = line[1:]  # remove ">"
                parts = header.split("$")
                if len(parts) != 4:
                    raise ValueError(f"Header not in expected format: {header}")
                accession, function, metabolic_function, organism = parts
                seq_lines = []
            else:
                seq_lines.append(line)
    
    # save last record
    if accession is not None:
        records.append({
            "accession": accession,
            "function": function,
            "metabolic_function": metabolic_function,
            "organism": organism,
            "aa_sequence": "".join(seq_lines)
        })

    return pd.DataFrame(records)


def filter_df_to_fasta(df: pd.DataFrame, path: str, function_value: str):
    """
    Filter DataFrame by function column and export filtered sequences to FASTA.
    """
    filtered_df = df[df["function"] == function_value]
    df_to_fasta(filtered_df, path)


def extract_header_info(fasta_file: str,
                        info: Literal["full", "accession", "dtype", "metabolic", "organism"] = "full",
                        dtype: Literal["as_set", "as_list"] = "as_set",
                    ) -> set[str] | list[str]:
    """
    Extract the header information from a FASTA file. 

    The header format is expected to be 
        "<accession>$<function>$<metabolic_function>$<organism>"
            each category, e.g. organism is separated underscorse instead of spaces

    For info, choose between: 
        - full: entire header 
        - accession: only the accession numbers
        - dtype: only the functions 
        - metabolic: only the metabolic functions
        - organism

    For dtype (data type), choose between:
        - as_set (default)
        - as_list 
    """

    if dtype == "as_set":
        infos = set()
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                header = record.description
                acc, func, meta, orga = header.split("$")
                if info == "full":
                    infos.add(header)
                elif info == "accession":
                    infos.add(acc)
                elif info == "function":
                    infos.add(func)
                elif info == "metabolic":
                    info.add(meta.replace("_", " "))
                else:
                    info.add(orga.replace("_", " "))
    else:
        infos = list()
        with open(fasta_file, "r") as f:
            for record in SeqIO.parse(f, "fasta"):
                header = record.description
                acc, func, meta, orga = header.split("$")
                if info == "full":
                    infos.append(header)
                elif info == "accession":
                    infos.append(acc)
                elif info == "function":
                    infos.append(func)
                elif info == "metabolic":
                    info.append(meta.replace("_", " "))
                else:
                    info.append(orga.replace("_", " "))

    return infos
                


# %%
def merge_fastas(input_fastas, output_fasta):
    """
    Merge FASTA files 
     -> remove duplicate sequences 
    Merge FASTA files 
     -> remove duplicate sequences 

    Parameters
    ----------
    input_fastas : list of str
        Paths to input FASTA files.
    output_fasta : str
        Path to merged FASTA file.
    """
    seen_ids = set()
    output_records = []

    for fasta in input_fastas:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id not in seen_ids:
                seen_ids.add(record.id)
                output_records.append(record)

    # write merged fasta
    seen_ids = set()
    output_records = []

    for fasta in input_fastas:
        for record in SeqIO.parse(fasta, "fasta"):
            if record.id not in seen_ids:
                seen_ids.add(record.id)
                output_records.append(record)

    # write merged fasta
    with open(output_fasta, "w") as out:
        SeqIO.write(output_records, out, "fasta")

    return output_fasta



def sanitize_fasta_headers(input_fasta, output_fasta):
    """
    Replace illegal characters in FASTA headers to avoid breaking MAFFT.
    Removes internal '>' and trims whitespace.
    """
    with open(output_fasta, "w") as out:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Remove ">" from anywhere in the header
            clean_header = record.description.replace(">", "_")
            
            # You may also want to replace spaces:
            clean_header = clean_header.replace(" ", "_")

            record.id = clean_header
            record.description = clean_header
            
            SeqIO.write(record, out, "fasta")


def remove_metabolic_function(input_fasta: str, output_fasta: str):
    """
    Removes the 3rd field from headers structured as:
        >acc$function$metabolic$organism
    Outputs:
        >acc$function$organism
    """
    with open(output_fasta, "w") as out:
        for rec in SeqIO.parse(input_fasta, "fasta"):
            header = rec.description

            # split on '$'
            parts = header.split("$")

            # if the header does not match expected format, keep it unchanged
            if len(parts) >= 4:
                # remove the 3rd part (index 2)
                new_parts = [parts[0], parts[1], parts[3]]
                new_header = "$".join(new_parts)
            else:
                # fallback: leave header unchanged
                new_header = header

            # write modified record
            out.write(f">{new_header}\n")
            seq = str(rec.seq)
            out.write("\n".join(seq[i:i+80] for i in range(0, len(seq), 80)))
            out.write("\n")





def clean_fasta_file(
    fasta_path: Path,
    species: str,
    tax_info: dict[str, str | int] | int,
    output_dir: Path,
    max_header_length: int = 80,
):
    """
    Clean FASTA headers, append taxid, rewrite FASTA in place,
    and write a JSON mapping old -> new headers.

    Output JSON:
        output_dir / "old_to_new_fasta_headers.json"
    """

    def clean_header(header: str) -> str:
        # 1) remove leading ">"
        h = header.lstrip(">")

        # 2) replace spaces with underscore
        h = header.replace(" ", "_")

        # 3) keep only safe characters
        # allowed: letters, numbers, _, ., -, |
        h = re.sub(r"[^A-Za-z0-9_.\-|]", "", h)

        # 4) trim underscores at ends
        h = h.strip("_")

        return h

    # taxid
    if isinstance(tax_info, int):
        taxid = tax_info
    else:
        if species not in tax_info:
            raise ValueError(f"Species '{species}' not found in species_to_taxid")

        taxid = tax_info[species]


    # read records
    records = list(SeqIO.parse(fasta_path, "fasta"))

    old_to_new: dict[str, str] = {}
    used_headers: set[str] = set()

    # process headers
    for record in records:
        old_header = record.description

        base = clean_header(old_header)

        # add taxid suffix
        new_header = f"{base}__{taxid}"

        # enforce max length
        if len(new_header) > max_header_length:
            cut = max_header_length - len(f"__{taxid}")
            base_cut = base[:cut]
            new_header = f"{base_cut}__{taxid}"

        # ensure uniqueness
        final_header = new_header
        counter = 1
        while final_header in used_headers:
            suffix = f"_dup{counter}"
            max_base_len = max_header_length - len(suffix)
            final_header = f"{new_header[:max_base_len]}{suffix}"
            counter += 1

        used_headers.add(final_header)
        old_to_new[old_header] = final_header

        # update record in place
        record.id = final_header
        record.name = final_header
        record.description = final_header

    # rewrite FASTA in place
    SeqIO.write(records, fasta_path, "fasta")


    # write JSON mapping
    output_dir.mkdir(parents=True, exist_ok=True)
    json_path = output_dir / "clean_fasta_headers.json"

    # merge with existing json if present
    if json_path.exists():
        with open(json_path, "r", encoding="utf-8") as f:
            try:
                existing = json.load(f)
            except json.JSONDecodeError:
                existing = {}
    else:
        existing = {}

    # add only new mappings
    for old, new in old_to_new.items():
        if old not in existing:
            existing[old] = new

    # write merged result
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(existing, f, indent=2)

    return old_to_new 

