from typing import Literal
import random
import warnings

def species_by_taxonomic_rank_picker(
    tax_ids: (list[int]) | (set[int]),
    by_rank: Literal["order", "family", "genus"] = "order",
    return_n_species: int = 50
) -> list[int]:
    """
    Randomly select species (represented by their taxonomic IDs) stratified by a 
    chosen taxonomic rank.

    Parameters
    ----------
    tax_ids : list[int]
        A list of NCBI taxonomic IDs representing species.

    by_rank : {"order", "family", "genus"}, optional
        Taxonomic rank used for grouping the species before sampling.
        Default is "order".

    return_n_species : int, optional
        Total number of species to return. Sampling is performed as follows:
        
        - Species are grouped by the chosen rank.
        - One species is randomly chosen from each rank (without replacement).
        - This process repeats until the requested number is reached or until all
          species are exhausted.
        
        Default: 50.

    Returns
    -------
    list[int]
        A list of selected taxonomic IDs.

    Raises
    ------
    ValueError
        If ``return_n_species`` is larger than the number of available species.

    Notes
    -----
    - The function assumes all IDs correspond to valid “species-level” taxonomic
      nodes in the NCBI taxonomy database.
    - Species that do not have the requested rank available will be ignored (with warning).
    """

    if return_n_species > len(tax_ids):
        raise ValueError("Total number of species to return must be smaller than number of parsed tax ids")
    
    tax_ids = set(tax_ids)

    # rank_species_dict = {}
    # for tid in tax_ids:
    #     rank = ncbi.

