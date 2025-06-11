#!/bin/bash

# Create a directory for the human genome and move into it
mkdir -p HumanGenome
cd HumanGenome

# Download genome sequence and annotation files for GRCh38 (Gencode v29)
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz

# Unzip the downloaded files and rename them for convenience
gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa Genome.fa

gunzip gencode.v29.annotation.gtf.gz
mv gencode.v29.annotation.gtf Genome.gtf

# Build STAR genome index
GENOMEDIR="$(pwd)"
mkdir -p "$GENOMEDIR/STAR"

STAR --runThreadN 23 \
     --runMode genomeGenerate \
     --genomeDir "$GENOMEDIR/STAR" \
     --genomeFastaFiles "$GENOMEDIR/Genome.fa" \
     --sjdbGTFfile "$GENOMEDIR/Genome.gtf"

# Go back to the parent directory
cd ..

# Prepare viral genome for BWA
mkdir -p ViralGenome
cd ViralGenome

# Download viral genome data
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
mv viral.1.1.genomic.fna viral_genome.fa

# Index the viral genome using BWA
bwa index -a bwtsw viral_genome.fa