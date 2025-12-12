from typing import Literal


def random_species_by_taxonomic_rank_picker(
        species: list[str | int], by_rank: Literal["order", "family", "genus"] = "order", 
        number_of_species: int= 50
        ) -> list[str]:
    """
    Take a list of species (either as scientific names or by taxonomic number) 
    and return Selects a specified number of species randomly from a given list. 
    """
    





