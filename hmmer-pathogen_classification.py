import os

def run_hmmer(query_faa, hmm_db, output_file):
    """Runs HMMER for protein family/domain search."""
    cmd = f"hmmsearch --tblout {output_file} {hmm_db} {query_faa}"
    os.system(cmd)

# Example Usage
run_hmmer("query_proteins.faa", "Pfam-A.hmm", "hmmer_results.tbl")

# Script written by Chukwuma Winner Obiora, Bioinformatician at University of Augsburg/Klinikum
