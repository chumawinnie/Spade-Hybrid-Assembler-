import os
import subprocess
import glob

# Define input directories
illumina_dir = "/mnt/share_uni_lehrstuhldaten/ftp/jandata/viral_dataset_from_Kevin_Dennehy_UKA/data_from_Kevin_26.8.2024"
nanopore_dir = os.path.expanduser("~/assembly_spades-project/trimmed_reads")
output_dir = os.path.expanduser("~/assembly_output")

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

# Get list of Illumina R1 fastq.gz files
illumina_samples = glob.glob(os.path.join(illumina_dir, "*_1.fastq.gz"))

# Process each sample
for illumina_r1 in illumina_samples:
    # Extract sample name (e.g., TS195, TS238, etc.)
    base_name = os.path.basename(illumina_r1)
    sample_parts = base_name.split("_")
    
    if len(sample_parts) < 2:
        print(f"⚠️ Skipping unexpected filename format: {base_name}")
        continue

    sample_name = sample_parts[1]

    # Construct paths
    read1_path = illumina_r1
    read2_path = read1_path.replace("_1.fastq.gz", "_2.fastq.gz")

    # Correct Nanopore naming pattern to match actual filenames
    nanopore_reads = os.path.join(nanopore_dir, f"merged_{sample_name}_nanopore_reads_trimmed.fastq.gz")

    # Output directories
    spades_out = os.path.join(output_dir, f"{sample_name}_spades")
    quast_out = os.path.join(output_dir, f"{sample_name}_quast")

    # Check if all required files exist
    missing_files = []
    if not os.path.exists(read1_path):
        missing_files.append(f"Illumina Read1: {read1_path}")
    if not os.path.exists(read2_path):
        missing_files.append(f"Illumina Read2: {read2_path}")
    if not os.path.exists(nanopore_reads):
        missing_files.append(f"Nanopore Reads: {nanopore_reads}")

    if missing_files:
        print(f"\n⚠️ Missing files for {sample_name}. Skipping...")
        for missing in missing_files:
            print(f"   🔴 MISSING: {missing}")
        continue  # Skip this sample

    # Debugging: Print verified paths
    print(f"\n✅ Processing sample: {sample_name}")
    print(f"   📂 Illumina Read1: {read1_path}")
    print(f"   📂 Illumina Read2: {read2_path}")
    print(f"   📂 Nanopore Reads: {nanopore_reads}")

    # Run SPAdes Hybrid Assembly using Conda environment
    spades_cmd = (
        f"conda run -n assembly spades.py "
        f"--pe1-1 {read1_path} --pe1-2 {read2_path} "
        f"--nanopore {nanopore_reads} -o {spades_out}"
    )

    print(f"\n🔄 Running SPAdes for {sample_name}...")
    subprocess.run(spades_cmd, shell=True, executable="/bin/bash", check=True)

    # Run QUAST on assembled genome using Conda environment
    quast_cmd = f"conda run -n quast_env quast.py {spades_out}/scaffolds.fasta -o {quast_out}"
    print(f"\n🔄 Running QUAST for {sample_name}...")
    subprocess.run(quast_cmd, shell=True, executable="/bin/bash", check=True)

    print(f"\n✅ Done processing {sample_name}!")

print("\n🎉 All samples processed!")

