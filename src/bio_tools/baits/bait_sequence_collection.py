#%%
import pandas as pd
import openpyxl
from typing import Dict, List

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
            organism = row["organism"]
            long_id = f"{accession}${function}${metabolic_function}${organism}"
            fasta_ids.append(long_id)

        df["fasta_id"] = fasta_ids
        return df

    def _get_mapping(self, key_col: str, value_col: str, ) -> Dict[str, List[str]]:
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

        mapping: Dict[str, List[str]] = {}
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

baits = BaitSequenceCollection(path=PATH_bait_seq_col)

COLORS_2ODDs = {
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
    "2ODD23" : "#fff2cc",
    "LFS" : "#f4e8c3",
    "2OG1" : "#ffe599",
    "C2'H" : "#ffd966",
    "F6'H" : "#dbc7b0",
    "S8H" : "#ffc466",
    "GSLOH" : "#e6bd5a",
    "GRS" : "#c8b478",
    "TIIAS" : "#c1b9a0",
    "D4H" : "#bb9e9e",
    "BX6" : "#e0bbbb",
    "FNSI" : "#f4cccc",
    "FNSI_F3H" : "#d6b3b3",
    "FNSI_FLS" : "#dcb8b8",
    "F3H" : "#f4cccc",
    "FLS_F3H" : "",
    "H6H" : "#e9d0db",
    "IDS" : "#e9d0db",
    "SLC" : "#cfe2f3",
    "GIM" : "#d0d2e5",
    "M2H" : "#c27ba0",
    "M2H_weak" : "#e1afbc",
    "DMR6" : "#c2d3e2",
    "S5H" : "#abbbc9",
    "S3H" : "#95a3af",
    "FLS" : "#b4a7d6",
    "LDOX" : "#8e7cc3",
    "JOX" : "#6fa8dc",
    "ACCO" : "#3d85c6",
    "T6OD" : "#3371a8",
    "COD" : "#316ca2",
    "SRG" : "#316a9f",
    "LBO" : "#2c6190"
}
# %%
