#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
})

# Function to merge genomic regions
merge_genomic_regions <- function(file_path, distance_threshold = 100) {
  lines <- readLines(file_path)
  df <- data.frame(text = lines, stringsAsFactors = FALSE) %>%
    separate(text, into = c("chr", "range"), sep = "\\s+") %>%
    separate(range, into = c("start", "end"), sep = "-")

  df$start <- as.integer(df$start)
  df$end <- as.integer(df$end)
  df <- df[order(df$chr, df$start), ]

  merge_intervals <- function(group) {
    merged <- list()
    current <- group[1, ]
    for (i in 2:nrow(group)) {
      row <- group[i, ]
      if ((row$start - current$end) <= distance_threshold) {
        current$end <- max(current$end, row$end)
      } else {
        merged <- append(merged, list(current))
        current <- row
      }
    }
    merged <- append(merged, list(current))
    do.call(rbind, merged)
  }

  merged_df <- df %>%
    group_by(chr) %>%
    group_modify(~ if (nrow(.x) > 1) merge_intervals(.x) else .x) %>%
    ungroup()

  merged_df$length <- merged_df$end - merged_df$start
  merged_df <- merged_df[order(-merged_df$length), ]
  return(as.data.frame(merged_df))
}

# Function to process one sample
process_sample <- function(sam) {
  paired_end <- file.exists(paste0("./fastq/", sam, "_2.fastq"))

  # STAR alignment
  star_cmd <- paste(
    "STAR --runThreadN 8 --runMode alignReads --genomeDir Ref_Genome/",
    "--sjdbGTFfile Ref_Genome/Genome.gtf --readFilesIn",
    paste0("./fastq/", sam, "_1.fastq"),
    if (paired_end) paste0(" ./fastq/", sam, "_2.fastq") else "",
    "--outFilterMismatchNmax 2 --outFileNamePrefix", sam,
    "--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate"
  )
  system(star_cmd)

  # Remove intermediate BAM
  system(paste("rm -f", paste0(sam, "Aligned.sortedByCoord.out.bam")))

  # BWA alignment
  bwa_cmd <- paste(
    "bwa mem -t 8 ViralGenome/viral_genome.fa",
    paste0(sam, "Unmapped.out.mate1"),
    if (paired_end) paste0(" ", sam, "Unmapped.out.mate2") else "",
    ">", paste0(sam, "_aln.sam")
  )
  system(bwa_cmd)

  # Convert SAM to BAM
  system(paste("samtools view -Sb -h", paste0(sam, "_aln.sam"), ">", paste0(sam, "_aln.bam")))

  # Count viral reads
  system(paste("samtools view", paste0(sam, "_aln.bam"),
               "| cut -f3 | sort | uniq -c | awk '{if ($1>=10000) print $0}' >",
               paste0(sam, "_viruses_count.txt")))

  # Sort BAM
  system(paste("samtools sort", paste0(sam, "_aln.bam"), ">", paste0(sam, "_aln_sorted.bam")))

  # Generate continuous regions
  depth_cmd <- paste(
    "samtools depth", paste0(sam, "_aln_sorted.bam"),
    "| awk '{if ($3>=2) print $0}'",
    "| awk '{ if ($2!=(ploc+1)) {if (ploc!=0){printf(\"%s %d-%d\\n\",$1,s,ploc);} s=$2} ploc=$2; }'",
    ">", paste0(sam, "_continuous_region.txt")
  )
  system(depth_cmd)

  # Merge regions
  merged_data <- merge_genomic_regions(paste0(sam, "_continuous_region.txt"))
  write.table(merged_data, paste0(sam, "_continuous_region_merged.txt"), quote = FALSE, row.names = FALSE)

  # Compute viral genome lengths
  system("awk '/^>/ {if (seqlen) {split(header, a, \" \"); gsub(\">\", \"\", a[1]); print a[1], seqlen;} header=$0; seqlen=0; next} {seqlen += length($0)} END {split(header, a, \" \"); gsub(\">\", \"\", a[1]); print a[1], seqlen;}' ViralGenome/viral_genome.fa > viral_lengths.txt")

  # Load viral genome lengths
  viral_lengths <- read.table("viral_lengths.txt", header = FALSE, sep = " ", col.names = c("chr", "size"))

  # Merge with lengths
  viral_merged_with_length <- merge(merged_data, viral_lengths, by = "chr")
  viral_merged_with_length$ratio <- 100 * (viral_merged_with_length$length / viral_merged_with_length$size)

  # Calculate cutoff based on FASTQ read count
  fastq_file <- paste0("./fastq/", sam, "_1.fastq")
  fastq_lines <- as.numeric(system(paste("wc -l <", shQuote(fastq_file)), intern = TRUE))
  cutoff <- (fastq_lines / 2) * 0.4

  # Load viral read count
  virus_data <- fread(paste0(sam, "_viruses_count.txt"), header = FALSE, skip = 1, col.names = c("reads", "chr"))

  # Merge with length info
  viral_merged_with_length_reads <- merge(viral_merged_with_length, virus_data, by = "chr")

  viral_merged_with_length_reads_ft=viral_merged_with_length_reads[viral_merged_with_length_reads$length <= viral_merged_with_length_reads$size,]


  # Apply filtering
  viral_merged_with_filtered_length <- viral_merged_with_length_reads_ft[
    viral_merged_with_length_reads_ft$ratio > 50 | viral_merged_with_length_reads_ft$reads > cutoff, ]

  # Save final list
  write.table(viral_merged_with_filtered_length, paste0(sam, "_final_virus_list.txt"), row.names = FALSE, quote = FALSE)

  # Clean up
  system(paste("rm -f", paste0(sam, "_ids.txt"), paste0(sam, "_continuous_region_500.txt"), paste0(sam, "_viruses_names.txt")))
}

# Main execution
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: script.R <sample_list.txt>")
}
sample_file <- args[1]
samples <- readLines(sample_file)

for (sam in samples) {
  message("Processing sample: ", sam)
  process_sample(sam)
}
