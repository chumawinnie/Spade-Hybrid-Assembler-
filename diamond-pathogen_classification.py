import os

def run_diamond(query_faa, db, output_file):
    """Runs DIAMOND BLASTp for fast protein classification."""
    cmd = f"diamond blastp -d {db} -q {query_faa} -o {output_file} --outfmt 6"
    os.system(cmd)

# Example Usage
run_diamond("query_proteins.faa", "nr", "diamond_results.txt")

# Script written by Chukwuma Winner Obiora, Bioinformatician at University of Augsburg/Klinikum
