import os
import subprocess
import glob

# Define directories
filtered_dir = os.path.expanduser("~/assembly_spades-project/filtered_reads")
minimap_output_dir = os.path.expanduser("~/assembly_spades-project/minimap_output")
medaka_output_dir = os.path.expanduser("~/assembly_spades-project/medaka_output")

# Ensure Medaka output directory exists
os.makedirs(medaka_output_dir, exist_ok=True)

# Get all FASTQ and FASTA files
filtered_fastq_files = glob.glob(os.path.join(filtered_dir, "*_nanopore_reads_filtered.fastq"))
filtered_fasta_files = glob.glob(os.path.join(filtered_dir, "*.fasta"))

# Ensure all FASTA index (.fai) and minimap2 index (.mmi) files are created
for fasta_file in filtered_fasta_files:
    fai_file = f"{fasta_file}.fai"
    mmi_file = f"{fasta_file}.map-ont.mmi"

    sample_name = os.path.basename(fasta_file).replace(".fasta", "")

    # Create .fai index if missing
    if not os.path.exists(fai_file):
        print(f"🔄 Creating FAI index for {sample_name}...")
        try:
            subprocess.run(f"samtools faidx {fasta_file}", shell=True, check=True)
            print(f"✅ FAI index created: {fai_file}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error creating FAI index for {sample_name}: {e}")
            continue  # Skip this sample

    # Create .mmi index if missing
    if not os.path.exists(mmi_file):
        print(f"🔄 Creating Minimap2 MMI index for {sample_name}...")
        try:
            subprocess.run(f"minimap2 -I 16G -x map-ont -d {mmi_file} {fasta_file}", shell=True, check=True)
            print(f"✅ Minimap2 index created: {mmi_file}")
        except subprocess.CalledProcessError as e:
            print(f"❌ Error creating MMI index for {sample_name}: {e}")
            continue  # Skip this sample

# Process each sample for Medaka
for fastq_file in filtered_fastq_files:
    sample_name = os.path.basename(fastq_file).replace("_nanopore_reads_filtered.fastq", "")

    fasta_file = os.path.join(filtered_dir, f"{sample_name}.fasta")
    output_dir = os.path.join(medaka_output_dir, f"{sample_name}_polished")

    # Ensure corresponding FASTA file exists
    if not os.path.exists(fasta_file):
        print(f"❌ ERROR: Missing FASTA file for {sample_name}. Skipping...")
        continue

    print(f"\n🔄 Running Medaka for {sample_name}...")

    # Run Medaka consensus
    medaka_cmd = f"medaka_consensus -i {fastq_file} -d {fasta_file} -o {output_dir} -m r1041_e82_400bps_sup_v5.0.0"
    try:
        subprocess.run(medaka_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ Medaka completed for {sample_name}. Polished assembly saved in: {output_dir}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running Medaka for {sample_name}: {e}")

print("\n🎉 All samples processed successfully with Medaka!")
