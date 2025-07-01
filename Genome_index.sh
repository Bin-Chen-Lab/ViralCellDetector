#!/bin/bash

set -euo pipefail

# Helper function for error messages
error_exit() {
  echo "[ERROR] $1" >&2
  exit 1
}

# Check for required tools
for cmd in wget STAR bwa; do
  command -v $cmd >/dev/null 2>&1 || error_exit "$cmd not found. Please install it before running the script."
done

# Input file check
INPUT_FILE="Genome_file.txt"
if [[ ! -f $INPUT_FILE ]]; then
  error_exit "Input file $INPUT_FILE not found."
fi

# Read values from the file
FASTA_URL=$(sed -n '1p' "$INPUT_FILE") || error_exit "Failed to read FASTA URL."
GTF_URL=$(sed -n '2p' "$INPUT_FILE") || error_exit "Failed to read GTF URL."

# Create working directory
mkdir -p "Ref_Genome"
cd "Ref_Genome" || error_exit "Failed to enter Ref_Genome directory."

# Download files
wget "$FASTA_URL" -O Genome.fa.gz || error_exit "Failed to download genome FASTA."
wget "$GTF_URL" -O Genome.gtf.gz || error_exit "Failed to download genome GTF."

# Unzip
gunzip Genome.fa.gz || error_exit "Failed to unzip Genome.fa.gz"
gunzip Genome.gtf.gz || error_exit "Failed to unzip Genome.gtf.gz"

# Build STAR index
STAR --runThreadN 23 --runMode genomeGenerate --genomeDir ./ --genomeFastaFiles Genome.fa --sjdbGTFfile Genome.gtf || error_exit "STAR indexing failed."

echo "[INFO] Reference genome indexing completed."

# Return to parent directory
cd .. || error_exit "Failed to return to parent directory."

# Prepare viral genome
mkdir -p ViralGenome
cd ViralGenome || error_exit "Failed to enter ViralGenome directory."

wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz || error_exit "Failed to download viral genome."
echo "[INFO] Viral genome downloaded."

gunzip viral.1.1.genomic.fna.gz || error_exit "Failed to unzip viral genome."
mv viral.1.1.genomic.fna viral_genome.fa || error_exit "Failed to rename viral genome."

bwa index -a bwtsw viral_genome.fa || error_exit "BWA indexing failed."

echo "[INFO] Viral genome indexing completed."