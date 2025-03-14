import os

def run_blast(query_fasta, db, output_file):
    """Runs BLASTn for nucleotide sequence homology search."""
    cmd = f"blastn -query {query_fasta} -db {db} -out {output_file} -outfmt 6"
    os.system(cmd)

# Example Usage
run_blast("assembled_genome.fasta", "nt", "blast_results.txt")

# Script written by Chukwuma Winner Obiora, Bioinformatician at University of Augsburg/Klinikum
