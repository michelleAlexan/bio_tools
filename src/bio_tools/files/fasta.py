import pandas as pd

def df_to_fasta(df: pd.DataFrame, path: str):
    """
    Export a DataFrame to FASTA format with header:
    >accession$function$metabolic_function$organism
    """
    with open(path, "w") as f:
        for _, row in df.iterrows():
            header_parts = [
                str(row["accession"]),
                str(row["function"]),
                str(row["metabolic_function"]),
                str(row["organism"]).replace(" ", "_")  # optional: replace spaces
            ]
            header = ">" + "$".join(header_parts)
            
            # wrap sequence at 80 chars
            seq = row["aa_sequence"]
            wrapped_seq = "\n".join(seq[i:i+80] for i in range(0, len(seq), 80))
            
            f.write(f"{header}\n{wrapped_seq}\n")

import pandas as pd

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