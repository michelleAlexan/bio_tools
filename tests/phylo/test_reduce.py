#%%
import pytest

from ete4 import PhyloTree
from pathlib import Path
import json
from Bio import SeqIO

from bio_tools.phylo.reduce import reduce_tree_to_family_representatives


# %%
expected_fam_red_mapping = {

    # Malvaceae
    "Gobar.D05G137100.1__3634": [
        "lcl_NW_019167982.1_cds_XP_022777201.1_32743__66656", 
        "Gorai.011G070500.1_PACid-26808224__29730",
        "Gobar.D10G070100.1__3634", 
        "Thecc.06G101600.1__3641",
        "OMO70940__210143", 
        "lcl_NW_019167982.1_cds_XP_022777202.1_32744__66656", 
        "Gobar.D05G137100.1__3634", 
        "Gorai.009G137900.1_PACid-26768904__29730"
    ], 



    # family Cucurbitaceae
    "Tcord05348.2__703396": [
        "Amacr27099.1__386130", 
        "Pparo03366.2__386216", 
        "Ssyna04510.1__386148", 
        "Tcord05348.2__703396",
        "Tmedo05871.1__3035263", 
        "Ascan05317.1__388286",
        "Ptrip05575.1__540787", 
        "Ptrip05574.1__540787",
        "Bdioi18627.1__3652",
        "Slaci04436.1__1132105",
        "Sdipt36255.1__1132104",
        "Schir26374.1__1132094",
        "Sunda20924.1__1134498", 
        "Spach05988.1__515875",
        "Sbade31434.1__386234",
        "Slanc31382.1__1132131",
        "Mfaba06305.1__329128",
        "Ccord05719.1__1568992",
        "CM029400.1.CM029400.1.g402.t1__3670",
        "Clana02077.1__3654",
        "Bsimp06597.1__388288",
        "Mjeff06939.1__515865",
        "Csagi25158.1__217630",
        "Cvari06775.1__869952",
        "CsatW800767.3__3659",
        "CsatW800766.3__3659",
        "Dsoco00756.3__229694"
    ], 

    # genus Fagopyrum
    "CM008280.1.CM008280.1.g987.t1__62330": [
        "TRINITY_DN135656_c1_g2_i1__516549", 
        "CM008280.1.CM008280.1.g987.t1__62330",
        "TRINITY_DN142313_c0_g1_i1__516549"
    ],

        # genus Quercus
    "Qurub.10G093800.1__3512": [
        "Qurub.10G093800.1__3512", 
        "Qurub.10G093900.1__3512", 
        "lcl_NC_044913.1_cds_XP_030937519.1_43670__97700"
    ]
}



# %%

def test_reduce_tree_to_family_representatives(tmp_path):
    t_path = Path(__file__).parents[1] / "data" / "test_reduce_to_fam_level.nwk"
    f_path = Path(__file__).parents[1] / "data" / "test_reduce_to_fam_level.fasta"

    t = PhyloTree(str(t_path), sp_naming_function=lambda name: name.split('__')[-1])
    t.annotate_ncbi_taxa(taxid_attr='species')


    result, path_filtered_fasta, path_filtered_json = reduce_tree_to_family_representatives(
        t,
        str(f_path),
        str(tmp_path)
    )

    assert sorted([sorted(v) for v in result.values()]) == sorted([sorted(v) for v in expected_fam_red_mapping.values()])
    assert sorted(list(result.keys())) == sorted(list(expected_fam_red_mapping.keys()))


    assert path_filtered_fasta.exists()
    assert path_filtered_json.exists()

    records_pre_filtered = list(SeqIO.parse(f_path, "fasta"))
    records_post_filtered = list(SeqIO.parse(path_filtered_fasta, "fasta"))

    discarded_seqs = []
    for rep, seq in result.items():
        for member in seq:
            if member != rep:
                discarded_seqs.append(member)

    assert set([r.id for r in records_post_filtered]) == set([r.id for r in records_pre_filtered if r.id not in discarded_seqs])
    

