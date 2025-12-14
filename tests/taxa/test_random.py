import pytest

from bio_tools.taxa.random import species_by_taxonomic_rank_picker
["Arabidopsis thaliana", "Zea mays", "Oryza sativa"], 
[3702, 4577, 4530],

@pytest.mark.parametrize(
    "tax_ids, rank, return_n_species, expected, error",  
    [
        (
            [3702, 4577, 4530],
            "order",
            4, 
            None,
            ValueError
        ), 
    ]
)
def test_species_by_taxonomic_rank_picker(tax_ids, rank, return_n_species, expected, error):

    if error:
        with pytest.raises(ValueError):
            species_by_taxonomic_rank_picker(
                tax_ids=tax_ids, 
                by_rank=rank, 
                return_n_species=return_n_species
                )
    else:        
        result = species_by_taxonomic_rank_picker(
            tax_ids=tax_ids, 
            by_rank=rank, 
            return_n_species=return_n_species
            )

        assert result == expected