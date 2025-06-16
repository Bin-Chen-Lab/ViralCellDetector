#!/bin/bash

# Read input from external file
INPUT_FILE="Genome_file.txt"

# Read values from the file
FASTA_URL=$(sed -n '1p' $INPUT_FILE)
GTF_URL=$(sed -n '2p' $INPUT_FILE)

# Create a working directory and move into it
mkdir -p "Ref_Genome"
cd "Ref_Genome"

# Download genome sequence and annotation files
wget "$FASTA_URL" -O Genome.fa.gz
wget "$GTF_URL" -O Genome.gtf.gz

# Unzip and rename for consistency
gunzip Genome.fa.gz
gunzip Genome.gtf.gz

# Build STAR genome index
#GENOMEDIR="$(pwd)"
#mkdir -p "STAR"
STAR --runThreadN 23  --runMode genomeGenerate --genomeDir ./  --genomeFastaFiles Genome.fa  --sjdbGTFfile Genome.gtf

echo "Reference Genome indexing completed."
# Go back to the parent directory
cd ..

# Prepare viral genome for BWA
mkdir -p ViralGenome
cd ViralGenome

# Download viral genome data
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
echo "Viral Genome downloaded."
gunzip viral.1.1.genomic.fna.gz
mv viral.1.1.genomic.fna viral_genome.fa
# Index the viral genome using BWA
bwa index -a bwtsw viral_genome.fa

echo "Viral Genome indexing completed."