import os
import subprocess
import glob

# Define directories
filtered_dir = os.path.expanduser("~/assembly_spades-project/filtered_reads")
minimap_output_dir = os.path.expanduser("~/assembly_spades-project/minimap_output")

# Ensure output directory exists
os.makedirs(minimap_output_dir, exist_ok=True)

# Get all filtered Nanopore read files
filtered_files = glob.glob(os.path.join(filtered_dir, "*.fastq"))

# Process each file for conversion, Minimap2 alignment, BAM sorting, and indexing
for filtered_file in filtered_files:
    base_name = os.path.basename(filtered_file)  # Extract filename
    sample_name = base_name.replace("_nanopore_reads_filtered.fastq", "")  # Get sample ID
    
    # Define output file paths
    fasta_file = os.path.join(filtered_dir, f"{sample_name}.fasta")
    sam_file = os.path.join(minimap_output_dir, f"{sample_name}.sam")
    sorted_bam_file = os.path.join(minimap_output_dir, f"{sample_name}.sorted.bam")

    # Convert FASTQ to FASTA
    print(f"\n🔄 Converting {filtered_file} to {fasta_file}...")
    seqtk_cmd = f"seqtk seq -a {filtered_file} > {fasta_file}"
    try:
        subprocess.run(seqtk_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ FASTQ converted to FASTA: {fasta_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error converting FASTQ to FASTA for {sample_name}: {e}")
        continue  # Skip to the next sample

    # Run Minimap2 alignment
    print(f"\n🔄 Running Minimap2 for {sample_name}...")
    minimap_cmd = f"minimap2 -ax map-ont {fasta_file} {fasta_file} > {sam_file}"
    try:
        subprocess.run(minimap_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ Minimap2 completed for {sample_name}. Alignment saved: {sam_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running Minimap2 for {sample_name}: {e}")
        continue  # Skip to the next sample

    # Convert SAM to sorted BAM
    print(f"🔄 Sorting and converting SAM to BAM for {sample_name}...")
    sort_cmd = f"samtools sort -o {sorted_bam_file} {sam_file}"
    try:
        subprocess.run(sort_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ Sorted BAM file created: {sorted_bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error sorting BAM for {sample_name}: {e}")
        continue

    # Index the BAM file
    print(f"🔄 Indexing BAM file for {sample_name}...")
    index_cmd = f"samtools index {sorted_bam_file}"
    try:
        subprocess.run(index_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ BAM index created for {sample_name}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error indexing BAM for {sample_name}: {e}")

print("\n🎉 All samples processed successfully with FASTA conversion, Minimap2, sorted BAM, and BAM indexing!")
