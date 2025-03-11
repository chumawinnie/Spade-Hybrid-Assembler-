import os
import subprocess
import glob

# Define input and output directories
trimmed_dir = os.path.expanduser("~/assembly_spades-project/trimmed_reads")
filtered_dir = os.path.expanduser("~/assembly_spades-project/filtered_reads")

# Ensure the output directory exists
os.makedirs(filtered_dir, exist_ok=True)

# Get all trimmed Nanopore read files
trimmed_files = glob.glob(os.path.join(trimmed_dir, "*.fastq.gz"))

# Process each file with Filtlong
for trimmed_file in trimmed_files:
    base_name = os.path.basename(trimmed_file)  # Extract filename
    sample_name = base_name.replace("_nanopore_reads_trimmed.fastq.gz", "")  # Get sample ID
    
    # Define output file name
    filtered_file = os.path.join(filtered_dir, f"{sample_name}_nanopore_reads_filtered.fastq")
    
    # Run Filtlong
    filtlong_cmd = f"filtlong --min_length 1000 {trimmed_file} > {filtered_file}"
    print(f"\n🔄 Running Filtlong for {sample_name}...")
    
    try:
        subprocess.run(filtlong_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ Filtlong completed for {sample_name}. Output saved: {filtered_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running Filtlong for {sample_name}: {e}")

print("\n🎉 All samples processed with Filtlong!")
