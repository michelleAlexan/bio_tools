import pytest
from bio_tools.files.fasta import write_clean_fasta_with_taxid


def test_write_clean_fasta_with_taxid_correct_output(tmp_path):

    fasta_content = """>sp|Q9Y6K9|Some weird protein (isoform 1) !!!
ATGC
>gene 2 with spaces & symbols###
ATGCGG
"""

    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_content)

    species = "Homo sapiens"
    species_to_taxid = {"Homo sapiens": 9606}

    mapping, _, _  = write_clean_fasta_with_taxid(
        input_fasta_path=fasta_path,
        output_dir=tmp_path,
        scientific_sp_name=species,
        tax_info=species_to_taxid,
        max_header_length=80,
    )

    clean_fasta_path = tmp_path / "clean_test.fasta"
    json_path = tmp_path / "clean_fasta_headers.json"

    assert clean_fasta_path.exists()
    assert json_path.exists()

    expected_fasta_content = """>sp|Q9Y6K9|Some_weird_protein_isoform_1__9606
ATGC
>gene_2_with_spaces__symbols__9606
ATGCGG
"""

    expected_json = {
        "sp|Q9Y6K9|Some weird protein (isoform 1) !!!": "sp|Q9Y6K9|Some_weird_protein_isoform_1__9606",
        "gene 2 with spaces & symbols###": "gene_2_with_spaces__symbols__9606",
    }

    with open(clean_fasta_path) as f:
        fasta = f.read()

    assert fasta == expected_fasta_content
    assert mapping == expected_json


def test_write_clean_fasta_with_taxid_species_not_found(tmp_path):

    fasta_content = """>sp|Q9Y6K9|Some weird protein (isoform 1) !!!
ATGC
>gene 2 with spaces & symbols###
ATGCGG
"""

    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_content)

    species = "Homo sapiens"
    species_to_taxid = {"Pan troglodytes": 9593}  # Different species

    with pytest.raises(ValueError):
        write_clean_fasta_with_taxid(
            input_fasta_path=fasta_path,
            output_dir=tmp_path,
            scientific_sp_name=species,
            tax_info=species_to_taxid,
            max_header_length=80,
        )   
    
