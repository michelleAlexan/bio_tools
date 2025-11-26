from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
from io import StringIO
import pandas as pd
import taxoniq 
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO, Entrez
from io import StringIO
import pandas as pd
import re
import time
import subprocess
from pathlib import Path

def fetch_full_protein(accession, email="your.email@here.com"):
    """Fetch full protein sequence given an accession."""
    Entrez.email = email
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        seq_data = handle.read()
        handle.close()

        # Extract the sequence (remove FASTA header)
        seq = "".join(seq_data.split("\n")[1:])
        return seq
    except Exception:
        return None


def extract_accession(hit_id):
    """
    Parse accession from hit_id.
    BLAST hit IDs come in MANY forms:
       'ref|XP_1234.1|' 
       'sp|Q9LS23.1|XYZ'  
       'XP_1234567'
       'gi|12345|ref|NP_12345.1|'
    """
    # Try the common cases
    m = re.search(r"[A-Z]{1,3}_\d+\.\d+", hit_id)
    if m:
        return m.group(0)

    m = re.search(r"[OPQ][0-9][A-Z0-9]{3}[0-9](\.\d+)?", hit_id)
    if m:
        return m.group(0)

    # Fallback: return the whole hit ID
    return hit_id


def blast_species_restricted_df(
        fasta_file: str,
        species_taxid: int,
        email: str,
        database: str = "nr",
        program: str = "blastp",
        hitlist_size: int = 100,
        delay: float = 0.4
    ):
    """
    Run BLAST restricted to species and return:
    - All HSP alignments
    - Full subject protein sequences (efetch)
    as a pandas DataFrame.
    """
    rows = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        print(f"Running BLAST for: {record.id}")

        tax_query = f"txid{species_taxid}[ORGN]"

        # Call BLAST
        handle = NCBIWWW.qblast(
            program=program,
            database=database,
            sequence=str(record.seq),
            entrez_query=tax_query,
            hitlist_size=hitlist_size,
            descriptions=hitlist_size,
            alignments=hitlist_size,
            format_type="XML"
        )

        blast_xml = handle.read()
        handle.close()

        blast_record = NCBIXML.read(StringIO(blast_xml))

        # Parse hits
        for alignment in blast_record.alignments:
            accession = extract_accession(alignment.hit_id)

            # Fetch full sequence (cached per accession)
            full_seq = fetch_full_protein(accession, email=email)
            time.sleep(delay)

            for hsp_i, hsp in enumerate(alignment.hsps, start=1):

                rows.append({
                    "query_id": record.id,
                    "query_length": len(record.seq),

                    "hit_id": alignment.hit_id,
                    "accession": accession,
                    "hit_def": alignment.hit_def,
                    "hit_length": alignment.length,

                    "hsp_num": hsp_i,
                    "bit_score": hsp.bits,
                    "evalue": hsp.expect,
                    "identity": hsp.identities,
                    "alignment_length": hsp.align_length,
                    "positives": hsp.positives,
                    "gaps": hsp.gaps,

                    "query_start": hsp.query_start,
                    "query_end": hsp.query_end,
                    "hit_start": hsp.sbjct_start,
                    "hit_end": hsp.sbjct_end,

                    "aligned_query_seq": hsp.query,
                    "aligned_hit_seq": hsp.sbjct,

                    "aa_sequence": full_seq, 

                    "organism": taxoniq.Taxon(species_taxid).scientific_name
                })

    df = pd.DataFrame(rows)
    return df



# TODO: while the command works from the CLI, this pyhton function does not work..why tho?
# # makeblastdb -in /Users/michellealexander/projects/bait_sequence_collection/data/2ODDs/2ODD_baits_filtered_short.fasta -dbtype prot -out 2ODD_ref_db

def make_blast_db(input_fasta, output_folder, db_name="blast_ref_db", db_type="prot"):
    """
    Create a BLAST database from a FASTA file, saving all files in a specified folder.

    Parameters
    ----------
    input_fasta : str or Path
        Path to the input FASTA file.
    output_folder : str or Path
        Folder where the BLAST DB files will be written. Created if it does not exist.
    db_name : str
        Prefix for the BLAST database files.
    db_type : str
        'prot' for protein, 'nucl' for nucleotide.

    Returns
    -------
    Path
        Path to the BLAST database prefix.
    """

    input_fasta = Path(input_fasta)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)  # make folder if missing

    db_prefix = db_name

    cmd = [
        "makeblastdb",
        "-in", str(input_fasta),
        "-dbtype", db_type,
        "-out", str(db_prefix),
        "-parse_seqids"
    ]

    print("Running command:", " ".join(cmd))
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=output_folder)

    if result.returncode != 0:
        print("makeblastdb failed!")
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)

    print(f"BLAST database created: {db_prefix}")
    return db_prefix

