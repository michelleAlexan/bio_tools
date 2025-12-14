from typing import Literal
import random
from collections import defaultdict

from bio_tools.taxa.taxonomy import get_taxonomic_ranks
from bio_tools.utils.constraint import ensure_no_duplicates

def species_picker_by_taxonomic_rank(
    tax_ids: (list[int]) | (set[int]),
    rank: Literal["order", "family", "genus"] = "order",
    block_size: int = 50, 
    seed: int | None = None
) -> list[set]:
    """
    Randomly select species (represented by their taxonomic IDs) stratified by a 
    chosen taxonomic rank. Return blocks of randomly selected species. 

    E.g., 
    Blocks of 2 species for species 
        "Arabidopsis thaliana", "Zea mays", "Oryza sativa" ([3702, 4577, 4530])
          with chosen rank "order" (here 1 Brassicales, and 2 Poales)
        may result in 
            [{3702, 4577}, {4530}]

    Parameters
    ----------
    tax_ids : list[int]
        A list of NCBI taxonomic IDs representing species.

    by_rank : {"order", "family", "genus"}
        Taxonomic rank used for grouping the species before sampling.
        Default is "order".

    block_size : int
        Total number of species to return. Sampling is performed as follows:
        
        - Species are grouped by the chosen rank.
        - One species is randomly chosen from each rank (without replacement).
        - This process repeats until the requested number is reached or until all
          species are exhausted.
        
        Default: 50.

    Returns
    -------
    list[set]
        A list where each value represents a block of randomly selected species stratified by rank.

    Raises
    ------
    ValueError
        If ``block_size`` is larger than the number of available species.

    Notes
    -----
    - The function assumes all IDs correspond to valid “species-level” taxonomic
      nodes in the NCBI taxonomy database.
    """
    # ----------- helper function --------------
    def stratified_random_picker(
    grouped: dict[int, set[int]],
    block_size: int,
    seed: int | None = None,
    )  -> list[set[int]]:
        """
        Randomly pick species in blocks, maximizing rank diversity per block.

        Parameters
        ----------
        grouped
            Mapping of rank_id -> set of species taxids
        block_size
            Number of species per block
        seed
            Optional random seed for reproducibility

        Returns
        -------
        List[Set[int]]
            List of blocks of selected species
        """
        if seed:
            random.seed(seed)

        # Make a mutable copy
        remaining = {k: set(v) for k, v in grouped.items()}
        blocks: list[set[int]] = []

        while remaining:
            block: set[int] = set()
            ranks = list(remaining.keys())
            random.shuffle(ranks)

            for rank in ranks:
                if len(block) >= block_size:
                    break

                species = random.choice(tuple(remaining[rank]))
                block.add(species)
                remaining[rank].remove(species)

                if not remaining[rank]:
                    del remaining[rank]

            blocks.append(block)

        return blocks
    

    # ----------- enforce constraints --------------
    if block_size > len(tax_ids):
        raise ValueError(f"Total number of species {len(tax_ids)} must be smaller than provided number of tax ids {block_size}")
    tax_ids = ensure_no_duplicates(tax_ids)

    # ----------- mapping all rank values to provided species -----------
    # e.g., when provided the species "Arabidopsis thaliana", "Zea mays", "Oryza sativa"
    # the orders Brassicales (3699) and Poales (38820) map to 
    # grouped = {3699: {3702}, 38820: {4530, 4577}}

    species_rank_dict = get_taxonomic_ranks(taxa=tax_ids, rank=rank, return_scientific_notations=False)
    grouped = defaultdict(set)
    for k, v in species_rank_dict.items():
        grouped[v].add(k)
    
    grouped = dict(grouped)
    blocks = stratified_random_picker(grouped, block_size=block_size, seed=seed)

    return blocks

    

    
    

