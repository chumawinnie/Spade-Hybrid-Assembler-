# **Hybrid Genome Assembly & Pathogen Classification Pipeline**

## **Overview**
This repository contains a comprehensive bioinformatics workflow that integrates **hybrid genome assembly, variant calling, phylogenetic analysis, and pathogen classification**. The pipeline is designed for **bacterial and viral genomic studies** using **Illumina and Nanopore sequencing data**.

---
## **Table of Contents**
1. [Introduction](#introduction)
2. [Hybrid Genome Assembly](#hybrid-genome-assembly)
3. [Variant Calling](#variant-calling)
4. [Phylogenetic Analysis](#phylogenetic-analysis)
5. [Pathogen Classification](#pathogen-classification)
6. [Software & Dependencies](#software--dependencies)
7. [Installation & Setup](#installation--setup)
8. [Running the Pipeline](#running-the-pipeline)
9. [Results Interpretation](#results-interpretation)

---
## **Introduction**
This pipeline enables **de novo hybrid genome assembly** using **short-read Illumina** and **long-read Nanopore sequencing data**. The assembled genomes undergo **variant calling** to detect genomic variations, **phylogenetic analysis** to infer evolutionary relationships, and **pathogen classification** using homology search tools.

---
## **Hybrid Genome Assembly**
**Objective:** Generate high-quality genome assemblies by combining **high-accuracy short reads** with **long Nanopore reads** for better genome continuity.

**Tools Used:**
- **SPAdes**: Hybrid assembly of Illumina and Nanopore reads
- **Filtlong**: Quality filtering for Nanopore reads
- **Minimap2 & Racon**: Read correction and polishing
- **Medaka**: Additional polishing for Nanopore reads
- **QUAST**: Assembly quality assessment

### **Command Example:**
```bash
spades.py --pe1-1 illumina_R1.fastq.gz --pe1-2 illumina_R2.fastq.gz \
  --nanopore nanopore_reads.fastq.gz -o spades_output \
  --threads 16 --memory 64 --cov-cutoff off --isolate
```

---
## **Variant Calling**
**Objective:** Identify genomic variants (SNPs, indels, and structural variants) from assembled genomes.

**Tools Used:**
- **bcftools**: Variant calling and filtering
- **GATK**: HaplotypeCaller for SNP and indel detection

### **Command Example:**
```bash
bcftools mpileup -Ou -f reference.fasta aligned_reads.bam | bcftools call -mv -Ob -o variants.bcf
```

---
## **Phylogenetic Analysis**
**Objective:** Construct phylogenetic trees to infer evolutionary relationships among isolates.

**Tools Used:**
- **MAFFT**: Multiple sequence alignment
- **IQ-TREE**: Maximum likelihood phylogenetic tree construction

### **Command Example:**
```bash
mafft --auto sequences.fasta > aligned_sequences.fasta
iqtree -s aligned_sequences.fasta -m GTR+G -bb 1000 -nt AUTO
```

---
## **Pathogen Classification**
**Objective:** Identify and classify pathogens using homology search tools.

**Tools Used:**
- **BLAST**: Nucleotide/protein sequence alignment
- **DIAMOND**: Fast protein classification
- **HMMER**: Domain search and functional annotation

### **BLAST Command Example:**
```bash
blastn -query assembled_genome.fasta -db nt -out blast_results.txt -outfmt 6
```

### **DIAMOND Command Example:**
```bash
diamond blastp -d nr -q query_proteins.faa -o diamond_results.txt --outfmt 6
```

### **HMMER Command Example:**
```bash
hmmsearch --tblout results.tbl Pfam-A.hmm query_proteins.faa
```

---
## **Software & Dependencies**
| Tool        | Version  | Purpose |
|-------------|---------|---------|
| SPAdes      | 3.15.0  | Hybrid genome assembly |
| Minimap2    | 2.26    | Read alignment |
| Medaka      | 2.0.1   | Nanopore read polishing |
| QUAST       | 5.2.0   | Assembly quality assessment |
| bcftools    | 1.12    | Variant calling |
| MAFFT       | 7.490   | Multiple sequence alignment |
| IQ-TREE     | 2.1.4   | Phylogenetic tree construction |
| BLAST       | Latest  | Nucleotide/protein homology search |
| DIAMOND     | Latest  | Fast protein classification |
| HMMER       | Latest  | Hidden Markov Model search |

---
## **Installation & Setup**
Clone the repository and set up dependencies:
```bash
git clone https://github.com/chumawinnie/Spade-Hybrid-Assembler.git
cd Spade-Hybrid-Assembler
```

Create and activate the conda environments:
```bash
conda create -n assembly_env python=3.8
conda activate assembly_env
conda install -c bioconda spades minimap2 racon medaka quast bcftools mafft iqtree blast diamond hmmer
```

---
## **Running the Pipeline**
Run the entire pipeline step by step:
```bash
bash run_hybrid_assembly.py
bash run_variant_calling.py
bash run_phylogenetics.py
bash run_pathogen_classification.py
```

---
## **Results Interpretation**
- **Genome Assembly:** Check assembly statistics using QUAST.
- **Variant Calling:** Review SNPs and indels in the VCF output.
- **Phylogenetic Analysis:** Examine tree structures to infer relationships.
- **Pathogen Classification:** Identify potential pathogens using BLAST, DIAMOND, or HMMER results.

---
## **Future Work**
- Improve hybrid assembly by testing alternative polishing methods.
- Enhance pathogen classification with deep learning models.
- Automate full pipeline using Snakemake or Nextflow.

For questions or contributions, feel free to open an issue or submit a pull request.

