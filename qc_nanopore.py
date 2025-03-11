import os

# Define absolute paths
input_dir = "/mnt/share_uni_lehrstuhldaten/ftp/jandata/viral_dataset_from_Kevin_Dennehy_UKA/data_from_Hellen_21.2.2025__20240611_0706_NanoporeSeq/fastq_pass/merged_fastq_files"
output_dir = os.path.expanduser("~/assembly_spades-project")
trimmed_dir = os.path.join(output_dir, "trimmed_reads")
cleaned_dir = os.path.join(output_dir, "cleaned_reads")
fastqc_dir = os.path.join(output_dir, "fastqc_reports")

# Ensure output directories exist
os.makedirs(trimmed_dir, exist_ok=True)
os.makedirs(cleaned_dir, exist_ok=True)
os.makedirs(fastqc_dir, exist_ok=True)

# Define quality control parameters
quality_threshold = 7  # Minimum quality score
headcrop = 50          # Trim first 50 bases

# List all FASTQ files in the input directory
fastq_files = [f for f in os.listdir(input_dir) if f.endswith(".fastq.gz")]

# Process each file
for fastq in fastq_files:
    input_file = os.path.join(input_dir, fastq)
    
    print(f"\nProcessing: {input_file}")

    # Generate output filenames
    trimmed_fastq = os.path.join(trimmed_dir, fastq.replace(".fastq.gz", "_trimmed.fastq.gz"))
    cleaned_fastq = os.path.join(cleaned_dir, fastq.replace(".fastq.gz", "_cleaned.fastq.gz"))

    # Step 1: Trim adapters using Porechop
    print("Running Porechop for adapter removal...")
    cmd_porechop = f"porechop -i {input_file} -o {trimmed_fastq}"
    os.system(cmd_porechop)

    # Step 2: Quality filtering with NanoFilt
    print("Running NanoFilt for quality filtering...")
    cmd_nanofilt = f"python3 -m NanoFilt -q {quality_threshold} --headcrop {headcrop} < {trimmed_fastq} | gzip > {cleaned_fastq}"
    os.system(cmd_nanofilt)

    # Step 3: Run FastQC on the cleaned files
    print("Running FastQC on cleaned reads...")
    cmd_fastqc = f"fastqc -o {fastqc_dir} {cleaned_fastq}"
    os.system(cmd_fastqc)

    print(f"Finished processing: {cleaned_fastq}\n")

print("All files processed successfully!")
