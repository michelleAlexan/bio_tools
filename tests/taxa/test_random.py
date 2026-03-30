import pytest

from bio_tools.taxa.random import species_picker_by_taxonomic_rank
["Arabidopsis thaliana", "Zea mays", "Oryza sativa"], 
[3702, 4577, 4530],

@pytest.mark.parametrize(
    "tax_ids, rank, block_size, expected, error",  
    [
        (
            [3702, 4577, 4530],
            "order",
            4, 
            None,
            ValueError
        ), 
        (
            [3702, 4577, 4530],
            "order",
            2, 
            [{3702, 4577}, {4530}],
            None
        ), 
    ]
)
def test_species_picker_by_taxonomic_rank(tax_ids, rank, block_size, expected, error):

    if error:
        with pytest.raises(ValueError):
            species_picker_by_taxonomic_rank(
                tax_ids=tax_ids, 
                rank=rank, 
                block_size=block_size, 
                seed=42
                )
    else:        
        result = species_picker_by_taxonomic_rank(
            tax_ids=tax_ids, 
            rank=rank, 
            block_size=block_size, 
            seed=42
            )

        assert result == expected


