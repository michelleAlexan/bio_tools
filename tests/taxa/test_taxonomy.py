import pytest
from pathlib import Path

from bio_tools.taxa.taxonomy import (
    scientific_notation_to_tax_id, 
    ensure_taxa_level_is_lower_than_rank_level, 
    get_taxonomic_ranks
)

@pytest.mark.parametrize(
        "species, expected, warning",
[           
    ( # single value (not a list)
        "Arabidopsis thaliana", 
        3702,
        None

    ),  
    ( # list with all spiecies as string
        ["Arabidopsis thaliana", "Zea mays", "Oryza sativa"], 
        [3702, 4577, 4530],
        None
    ),  
    ( # a misspelled species name that is not provided in the up_to_date_scientific_notations_yaml
        ["Arabidopsis thalian", "Zea mays", "Oryza sativa"], 
        [4577, 4530],
        UserWarning
    ),  
        ( # a misspelled species name that is provided in the up_to_date_scientific_notations_yaml
        "Oryza sativ", 
        4530, 
        None
    ),  

],
)
def test_scientific_name_to_tax_id(species, expected, warning):
    yaml_path = Path(__file__).parents[1] / "data/example.yaml"

    # expected is an UserWarning (species not found) and is removed from returned list
    if warning:
        with pytest.warns(UserWarning):
            result = scientific_notation_to_tax_id(
                species=species,
                update_ncbi_db=False,
                up_to_date_scientific_notations_yaml=yaml_path,
            )
        #  check that invalid species was discarded 
        assert result == expected
    else:
        # expected is a normal value → compare outputs
        result = scientific_notation_to_tax_id(
            species=species,
            update_ncbi_db=False,
            up_to_date_scientific_notations_yaml=yaml_path,
        )
        assert result == expected



@pytest.mark.parametrize(
    "taxa, rank, val_error, key_error",
    [
        # invalid argument passed for 'rank' parameter
        (
            [1, 2, 3],
            "invalid_rank",
            ValueError,
            None
        ),
        # taxa contain species, genus, family and order, the provided rank is class → no error
        (
            [3702, 3701, 3700,3699],
            "class",
            None,
            None
        ),
        # taxa contains class & the provided rank is class → value error
        (
            [3398],
            "class",
            ValueError,
            None
        ),
    
        # taxa contains clade → key error
        (
            [3193],
            "class",
            None,
            KeyError,
        ),
    ],
)
def test_ensure_taxa_level_is_lower_than_rank_level(
    taxa, rank, val_error, key_error
):
    if val_error:
        with pytest.raises(ValueError):
            ensure_taxa_level_is_lower_than_rank_level(
                taxa=taxa, rank=rank
            )
    elif key_error:
        with pytest.raises(KeyError):
            ensure_taxa_level_is_lower_than_rank_level(
                taxa=taxa, rank=rank
            )
    else:
        ensure_taxa_level_is_lower_than_rank_level(taxa, rank)


@pytest.mark.parametrize(
    "taxa, rank, return_scientific_notation, expected, error",  
    [
        (
            "Arabidopsis", 
            "order",
            False,
            {
                3701: 3699
            },
            None,
        ),
        ( 
            "Arabidopsis", 
            "order",
            True,
            {
                'Arabidopsis': 'Brassicales'
            },
            None,
        ),

        ( 
            ["Arabidopsis thaliana", "Zea mays", "Oryza sativa"], 
            "genus",
            True,
            {
                "Arabidopsis thaliana": 'Arabidopsis', 
                "Zea mays": 'Zea', 
                "Oryza sativa": "Oryza", 
            },
            None,
        ),
    ]
)
def test_get_taxonomic_ranks(taxa, rank, return_scientific_notation, expected, error):
    if error:
        with pytest.raises(ValueError):
            get_taxonomic_ranks(
                taxa=taxa, 
                rank=rank, 
                return_scientific_notations=return_scientific_notation
            )
    else:
        result = get_taxonomic_ranks(
            taxa=taxa, 
                rank=rank, 
                return_scientific_notations=return_scientific_notation
            )
        assert result == expected

