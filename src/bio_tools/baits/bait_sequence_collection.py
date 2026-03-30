#%%
import pandas as pd
import openpyxl


# TODO: currently, I only have my characterized 2ODDs in the collection
# restrucutre class so that ingroup sequences that are not characterized 
# are well integrated!

class BaitSequenceCollection():
    """
    Represents a collection of ingroup and outgroup bait sequences for automated functional annotation.

    Expected Excel file:
      - sheet "ingroup"
      - sheet "outgroup"

    Each sheet should be a table with at least these columns (case-insensitive):
      - accession
      - function  
      - metabolic_function
      - organism

    The class builds:
      - self.ingroup_df
      - self.outgroup_df
      - self.combined_df   (ingroup + outgroup)
      - self.funcs         (set of all functions)
      - helper mappings: funcs -> fasta_ids, metabolic_function -> functions    """

    REQUIRED_COLS = {"accession", "function", "metabolic_function", "organism", "aa_sequence"}

    def __init__(self, path: str):
        """
        Parameters
        ----------
        path : str
            Path to the Excel workbook

        """

        try:
            self.ingroup_df = pd.read_excel(path, sheet_name="ingroup")
        except Exception as e:
            raise ValueError(f"Failed to read ingroup sheet: {e}")

        try:
            self.outgroup_df = pd.read_excel(path, sheet_name="outgroup")
        except Exception as e:
            raise ValueError(f"Failed to read outgroup sheet: {e}")

        # Validate required columns exist in both sheets
        missing_cols = self.REQUIRED_COLS - set(self.ingroup_df.columns)
        if missing_cols:
            raise ValueError(f"Missing required columns in ingroup sheet: {(missing_cols)}")

        # Add fasta_id (long format) 
        self.ingroup_df = self._add_fasta_id_long(self.ingroup_df)
        # self.outgroup_df = self._add_fasta_id_long(self.outgroup_df)

        # Set of all 2ODD functions
        self.funcs = set(self.ingroup_df["function"])

        # all unique fasta header ids
        self.fasta_ids = set(self.ingroup_df["fasta_id"])

        # Mappings
        self.funcs_fastaid_dict = self._get_mapping("function", "fasta_id")
        self.metabolic_funcs_dict = self._get_mapping("metabolic_function", "function")

    def _add_fasta_id_long(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Add a 'fasta_id' column using the format:
            <accession>$<function>$<metabolic_function>$<organism>

        If accession is missing/NaN/empty it will be replaced by a stable NONE{n} id per-dataframe,
        optionally prefixed by 'IN' or 'OUT' to keep ingroup/outgroup distinct.
        """
        df = df.copy()
        fasta_ids = []
        for _, row in df.iterrows():
            accession = row["accession"]
            function = row["function"]
            metabolic_function = row["metabolic_function"]
            organism = row["organism"].replace(" ", "_")
            long_id = f"{accession}${function}${metabolic_function}${organism}"
            fasta_ids.append(long_id)

        df["fasta_id"] = fasta_ids
        return df

    def _get_mapping(self, key_col: str, value_col: str, ) -> dict[str, list[str]]:
        """
        Return a mapping: key -> list(values) from self.combined_df.

        Parameters
        ----------
        key_col : str
            Column name to use as keys 
        value_col : str
            Column name to collect values for each key

        Returns
        -------
        dict
            mapping key -> list of values
        """

        mapping: dict[str, list[str]] = {}
        for i , row in self.ingroup_df[[key_col, value_col]].iterrows():
            key = row[key_col]
            val = row[value_col]
            if pd.isna(key) or str(key).strip() == "":
                print(f"WARNING: missing {key_col} at index {i}")

            if key not in mapping:
                mapping[key] = set()
            mapping[key].add(val)

        return mapping

    

#%%

PATH_bait_seq_col = "/Users/michellealexander/projects/bio_tools/src/bio_tools/baits/baits_copy.xlsx"
PATH_fasta = "/Users/michellealexander/projects/bait_sequence_collection/data/2ODDs/2ODD_baits_filtered.fasta"


char_2ODD_baits = BaitSequenceCollection(path=PATH_bait_seq_col)
from bio_tools.files.fasta import extract_header_info
fasta_headers = extract_header_info(PATH_fasta)

print(fasta_headers - char_2ODD_baits.fasta_ids)
print(char_2ODD_baits.fasta_ids - fasta_headers)

# %%

