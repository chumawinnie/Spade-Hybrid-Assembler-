import os
import subprocess
import glob

# Define directories
filtered_dir = os.path.expanduser("~/assembly_spades-project/filtered_reads")
racon_dir = os.path.expanduser("~/assembly_spades-project/racon_corrected")

# Ensure output directory exists
os.makedirs(racon_dir, exist_ok=True)

# Get all filtered Nanopore read files
filtered_files = glob.glob(os.path.join(filtered_dir, "*.fastq"))

# Process each file for correction
for filtered_file in filtered_files:
    base_name = os.path.basename(filtered_file)  # Extract filename
    sample_name = base_name.replace("_nanopore_reads_filtered.fastq", "")  # Get sample ID
    
    # Define output file paths
    sam_file = os.path.join(racon_dir, f"{sample_name}.sam")  # Minimap2 output
    racon_corrected_file = os.path.join(racon_dir, f"{sample_name}_racon_corrected.fastq")  # Racon output

    # Run Minimap2 alignment
    minimap_cmd = f"minimap2 -ax map-ont {filtered_file} {filtered_file} > {sam_file}"
    print(f"\n🔄 Running Minimap2 for {sample_name}...")
    
    try:
        subprocess.run(minimap_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ Minimap2 completed for {sample_name}. Alignment saved: {sam_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running Minimap2 for {sample_name}: {e}")
        continue  # Skip to next sample

    # Run Racon correction
    racon_cmd = f"racon {filtered_file} {sam_file} {filtered_file} > {racon_corrected_file}"
    print(f"\n🔄 Running Racon for {sample_name}...")
    
    try:
        subprocess.run(racon_cmd, shell=True, executable="/bin/bash", check=True)
        print(f"✅ Racon completed for {sample_name}. Corrected reads saved: {racon_corrected_file}")
    except subprocess.CalledProcessError as e:
        print(f"❌ Error running Racon for {sample_name}: {e}")

print("\n🎉 All samples processed with Racon!")
