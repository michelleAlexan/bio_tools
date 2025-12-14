#%%
from ete3 import NCBITaxa
from pathlib import Path
from collections.abc import Iterable
from typing import Literal
import warnings

from bio_tools.utils.io import load_yaml
from bio_tools.utils.constraint import internet_on, ensure_list, ensure_no_duplicates

ncbi = NCBITaxa()

RANK_ORDER = {
    "species": 0,
    "subspecies": 1,
    "genus": 2,
    "family": 3,
    "order": 4,
    "class": 5,
    "phylum": 6,
    "kingdom": 7,
    "superkingdom": 8,
}

UP_TO_DATE_SCIENTIFIC_NOTATIONS_YAML = Path("/Users/michellealexander/projects/bio_tools/config/up_to_date_species_name.yaml")
#%% 



def map_species_to_correct_names(scientific_names: (list[str] )| (str), 
                                 up_to_date_scientific_notations_yaml: Path):
    """
    Take a scientific name(s) that are misspelled 
        and take a path to a yaml file in which wrong names map to up-to-date (correct) names. 
    Return the same list with mapped values. 
    """
    up_to_date_scientific_notations_dict = load_yaml(up_to_date_scientific_notations_yaml)

    scientific_names = ensure_list(scientific_names)
    corrected_list = [up_to_date_scientific_notations_dict[s] if s in up_to_date_scientific_notations_dict.keys() else s for s in scientific_names]
    return corrected_list


def scientific_notation_to_tax_id(species:(list[str] )| (str), 
                              update_ncbi_db: bool = False, 
                              up_to_date_scientific_notations_yaml: Path | None = None
                            ) -> (list[int] | int):
    """
    Take a list of species, and return as list with corresponding taxon ids.
    
    For this, convert scientific name to taxonomic id using the ete3.NCBITaxa package. 
        If used for first time, a stable internet is required to fetch the taxonomy database.
            - > It will be saved in under '~/.etetoolkit/taxa.sqlite'
        If you haven't updated it for a while and want to fetch the up-to-date database, 
            - > set 'update_ncbi_db' parameter to True (may take 2 minutes)

    """

    def ete3_ncbi_approach(taxa):
        """
        Helper function. Fetch the species tax ids using the ete3.NCBITaxa approach. 
        """
        tax_ids = []
        for s in taxa:
            try:
                taxid = ncbi.get_name_translator([s])[s][0]
                tax_ids.append(taxid)
            except KeyError:
                warnings.warn(f"The species '{s}' was not found by ete3.NCBITaxa.get_name_translator and is thus discarded."
                            "Check for spelling mistakes or if species name is depricated.", UserWarning)
                continue

        return tax_ids
    
    # ----------- enforce constraints --------------
    # normalize to list
    species = ensure_list(species)
    print(species)
    # remove duplicates
    species = ensure_no_duplicates(species)
    print(species)


    # ----------- handle potential typos / depricated scientific names --------------
    # In case the list of scientific species name includes any typos or depricated species names, 
    # it is possible to manually create a yaml file containing the mapping:
    #    depricated_scientific_notation -> up_to_date_scientific_notation.
    if up_to_date_scientific_notations_yaml:
        species = map_species_to_correct_names(species, up_to_date_scientific_notations_yaml)


    # ---------- fetch tax ids using ete3.NCBITaxa approach --------------
    # this fetches the most up to date species - tax id mapping.
    if Path('~/.etetoolkit/taxa.sqlite').expanduser().exists():
        # update taxonomy database 
        if internet_on() and update_ncbi_db:
            ncbi.update_taxonomy_database()
        else:
            tax_ids = ete3_ncbi_approach(species)

    # todo: implement taxoniq offline method

    # --------- return -----------
    if len(tax_ids) == 1:
        return tax_ids[0]
    else:
        return tax_ids


def ensure_taxa_level_is_lower_than_rank_level(
    taxa: list[int],
    rank: Literal["class", "order", "family", "genus"],
) -> None:
    """
    Ensure that all taxa are strictly lower than the specified rank.

    Raises
    ------
    ValueError
        If any taxon has rank > specified rank
    """
    allowed_ranks = ["genus", "family", "order", "class"]

    if rank not in allowed_ranks:
        raise ValueError(
            f"Invalid rank '{rank}'. "
            f"Must be one of {allowed_ranks}."
        )
    
    target_rank_value = RANK_ORDER[rank]
    taxid_to_rank = ncbi.get_rank(taxa)

    for taxid, input_taxa_rank in taxid_to_rank.items():
        if input_taxa_rank not in RANK_ORDER:
            raise KeyError(
                f"TaxID {taxid} has unsupported rank '{input_taxa_rank}'. "
                f"Allowed ranks: {allowed_ranks}."
            )
        
        if RANK_ORDER[input_taxa_rank] >= target_rank_value:
            raise ValueError(
                f"TaxID {taxid} has rank '{input_taxa_rank}', "
                f"which is not lower than requested rank '{rank}'."
            )



def get_taxonomic_ranks(
    taxa: (list[int]) | (set[int]) | (list[str]) | (set[str]), 
    rank: Literal["genus", "family", "order", "class"], 
    return_scientific_notations: bool = False
) -> (dict[int, int]) | (dict[str, str]):
    """
    Fetch the specified taxonomic rank for each tax ID.

    Parameters
    ----------
    tax_ids : List[int]
        List of NCBI taxonomic IDs.
    rank : str
        Taxonomic rank to retrieve (e.g., "order", "family").

    Returns
    -------
    Dict[int, Optional[str]]
        A dictionary mapping:
        tax id or scientific notation of speciefied taxa -> tax id or scientific name of the specified rank.
        If the rank is not found for a tax_id, the value will be None.

    Notes
    -----
    - Uses ete3.NCBITaxa.get_lineage and get_rank to traverse the taxonomy tree.
    - The function warns for tax IDs not found in NCBI database.
    """


    # --------- convert scientific notation to tax ids ------------
    if isinstance(taxa, str) or isinstance(taxa, Iterable):
        if isinstance(taxa, Iterable) and all(isinstance(x, str) for x in taxa):
            taxa = scientific_notation_to_tax_id(taxa)


    # ----------- enforce constraints --------------
    taxa = ensure_list(taxa)
    taxa = ensure_no_duplicates(taxa)
    ensure_taxa_level_is_lower_than_rank_level(taxa, rank)

    # ---------- mapping all input values to corresponding ranks -------------
    
    result = {}
    for tid in taxa:
        # get full lineage 
        lineage = ncbi.get_lineage(tid)
        # get rank for all nodes in lineage
        ranks = ncbi.get_rank(lineage)  # dict {taxid: rank_name}
        # find the taxid that matches requested rank
        rank_taxid = next((taxid for taxid in lineage if ranks.get(taxid) == rank), None)
        if rank_taxid:
            # fetch the scientific notation, if requested
            if return_scientific_notations:
                tid_sn = next(iter(ncbi.get_taxid_translator([tid]).values()))
                rank_sn = next(iter(ncbi.get_taxid_translator([rank_taxid]).values()))
                result[tid_sn] = rank_sn
            else:
                result[tid] = rank_taxid

        else:
            warnings.warn(f"Tax ID {tid} does not have a rank '{rank}'", UserWarning)
            result[tid] = None

    return result


#%%


# TODO: TEST
def get_taxa_topology(taxa_ids: list[int]):
    """
    
    """
    ncbi = NCBITaxa()

    tree = ncbi.get_topology(taxa_ids, annotate=True)

    for node in tree.traverse():
        # Remove problematic characters from names
        if hasattr(node, "sci_name") and node.sci_name:
            node.name = node.sci_name.replace(":", "_").replace("(", "").replace(")", "").replace(" ", "_")
        elif not node.name:
            node.name = "NA"
    # Write a clean Newick file
    tree.write(outfile="plant_tree_clean.nwk", format=5)

# %%
