import pytest
from bio_tools.files.fasta import clean_fasta_file


def test_clean_fasta_file(tmp_path):

    fasta_content = """>sp|Q9Y6K9|Some weird protein (isoform 1) !!!
ATGC
>gene 2 with spaces & symbols###
ATGCGG
"""
    fasta_path = tmp_path / "test.fasta"
    fasta_path.write_text(fasta_content)

    species = "Homo sapiens"
    species_to_taxid = {"Homo sapiens": 9606}

    json_path = tmp_path / "clean_fasta_headers.json"
    mapping = clean_fasta_file(
        fasta_path=fasta_path,
        species=species,
        tax_info=species_to_taxid,
        output_dir=tmp_path,
        max_header_length=80,
    )

    expected_fasta_content = """>sp|Q9Y6K9|Some_weird_protein_isoform_1__9606
ATGC
>gene_2_with_spaces__symbols__9606
ATGCGG
"""    
    expected_fasta_path = tmp_path / "expected.fasta"
    expected_fasta_path.write_text(expected_fasta_content)

    assert json_path.exists()

    expected_json = {
        "sp|Q9Y6K9|Some weird protein (isoform 1) !!!" : "sp|Q9Y6K9|Some_weird_protein_isoform_1__9606",
        "gene 2 with spaces & symbols###": "gene_2_with_spaces__symbols__9606"
    }

    print(mapping)
    assert mapping == expected_json

    with open(fasta_path, "r", encoding="utf-8") as f:
        fasta = f.read()

    with open(expected_fasta_path, "r", encoding="utf-8") as f:
        expected = f.read()

    assert fasta == expected


