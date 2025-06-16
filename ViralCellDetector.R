#!/usr/bin/env Rscript

# Load required libraries
suppressMessages({
  require(data.table)
})

# Function to process one sample
process_sample <- function(sam) {
  paired_end <- file.exists(paste0("./fastq/", sam, "_2.fastq"))

  # Set up STAR alignment
  star_cmd <- paste(
    "STAR --runThreadN 8 --runMode alignReads --genomeDir Ref_Genome/",
    "--sjdbGTFfile Ref_Genome/Genome.gtf --readFilesIn",
    paste0("./fastq/", sam, "_1.fastq"),
    if (paired_end) paste0("./fastq/", sam, "_2.fastq") else "",
    "--outFilterMismatchNmax 2 --outFileNamePrefix", sam,
    "--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate"
  )
  system(star_cmd)

  # Remove unnecessary STAR BAM file
  system(paste("rm -f", paste0(sam, "Aligned.sortedByCoord.out.bam")))

  # Run BWA on unmapped reads
  bwa_cmd <- paste(
    "bwa mem -t 8 ViralGenome/viral_genome.fa",
    paste0(sam, "Unmapped.out.mate1"),
    if (paired_end) paste0(sam, "Unmapped.out.mate2") else "",
    ">", paste0(sam, "_aln.sam")
  )
  system(bwa_cmd)

  # Convert SAM to BAM
  system(paste("samtools view -Sb -h", paste0(sam, "_aln.sam"), ">", paste0(sam, "_aln.bam")))

  # Count unique viral alignments
  system(paste("samtools view", paste0(sam, "_aln.bam"), "| cut -f3 | sort | uniq -c | awk '{if ($1>=10000) print $0}' >", paste0(sam, "_viruses_count.txt")))

  # Sort BAM file
  system(paste("samtools sort", paste0(sam, "_aln.bam"), ">", paste0(sam, "_aln_sorted.bam")))

  # Get depth and continuous regions
  depth_cmd <- paste(
    "samtools depth", paste0(sam, "_aln_sorted.bam"),
    "| awk '{if ($3>=5) print $0}'",
    "| awk '{ if ($2!=(ploc+1)) {if (ploc!=0){printf(\"%s %d-%d\\n\",$1,s,ploc);} s=$2} ploc=$2; }'",
    ">", paste0(sam, "_continuous_region.txt")
  )
  system(depth_cmd)

  # Generate region ID info
  system(paste("cut -d ' ' -f2", paste0(sam, "_continuous_region.txt"),
               "| awk -F '[-]' '{print$1\"\\t\"$2\"\\t\"$2-$1}' >", paste0(sam, "_ids.txt")))

  # Filter long continuous regions
  system(paste("paste", paste0(sam,"_continuous_region.txt"), paste0(sam,"_ids.txt"),
               "|awk '{if($5 >=100)print$1\"\\t\"$3\"\\t\"$4\"\\t\"$5}' >", paste0(sam,"_continuous_region_100.txt")))

  # Extract virus names
  system(paste("awk '{print $2}'", paste0(sam, "_viruses_count.txt"), ">", paste0(sam, "_viruses_names.txt")))

  # Final virus list
  system(paste("fgrep -wf", paste0(sam, "_viruses_names.txt"), paste0(sam,"_continuous_region_100.txt"),
               ">", paste0(sam,"_final_virus_list.txt")))

  # Cleanup
  system(paste("rm -f", paste0(sam,"_ids.txt"), paste0(sam,"_continuous_region_100.txt"), paste0(sam, "_viruses_names.txt")))
}

# Main script execution
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
