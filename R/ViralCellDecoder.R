#!/usr/bin/env  Rscript

require(data.table)
require(plyr)
require(dplyr)
require(foreach)
require(parallel)
require(doParallel)
# Read in sample names from a file
filename <- commandArgs(trailingOnly=TRUE)
samples <- readLines(filename)

setwd("./")### set path as per you convenient

for (sam in samples) {
  #assuming all files are available in the same folder where the script is running
  #else add absolute path to following lines
  if (file.exists(paste0("./fastq/", sam, "_2.fastq"))) {
    system(paste("STAR --runThreadN 8 --runMode alignReads --genomeDir",
                 "HumanGenome/ --sjdbGTFfile",
                 "HumanGenome/Genome.gtf --readFilesIn",
                 paste0("./fastq/", sam, "_1.fastq"), paste0("./fastq/", sam, "_2.fastq"),
                 "--outFilterMismatchNmax 2 --outFileNamePrefix", sam,
                 "--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate"))

    system(paste("rm", paste0(sam, "Aligned.sortedByCoord.out.bam")))

    system(paste("bwa mem -t 8", "ViralGenome/viral_genome.fa",
                 paste0(sam, "Unmapped.out.mate1"), paste0(sam, "Unmapped.out.mate2"), ">", paste0(sam, "_aln.sam")))

    system(paste("samtools view -Sb -h", paste0(sam, "_aln.sam"), ">", paste0(sam, "_aln.bam")))

    system(paste("samtools view", paste0(sam, "_aln.bam"), "| cut -f3 | sort | uniq -c | awk '{if ($1>=10000) print $0}' >", paste0(sam, "_viruses_count.txt")))

    system(paste("samtools sort", paste0(sam, "_aln.bam"), ">", paste0(sam, "_aln_sorted.bam")))

    system(paste("samtools depth", paste0(sam, "_aln_sorted.bam"), "| awk '{if ($3>=5) print $0}' | awk '{ if ($2!=(ploc+1)) {if (ploc!=0){printf(\"%s %d-%d\\n\",$1,s,ploc);}s=$2} ploc=$2; }' >", paste0(sam, "_continuous_region.txt")))

    system(paste("cut -d ' ' -f2", paste0(sam, "_continuous_region.txt"), "| awk -F '[-]' '{print$1\"\\t\"$2\"\\t\"$2-$1}' >", paste0(sam, "_ids.txt")))

    system(paste("paste", paste0(sam,"_continuous_region.txt"), paste0(sam,"_ids.txt"),"|awk '{if($5 >=100)print$1\"\\t\"$3\"\\t\"$4\"\\t\"$5}' >", paste0(sam,"_continuous_region_100.txt")))

    system(paste("awk '{print $2}'", paste0(sam, "_viruses_count.txt"), ">",   paste0(sam, "_viruses_names.txt")))

    system(paste("fgrep -wf", paste0(sam, "_viruses_names.txt"), paste0(sam,"_continuous_region_100.txt"), ">", paste0(sam,"_final_virus_list.txt")))

    system(paste("rm", paste0(sam,"_ids.txt"), paste0(sam,"_continuous_region_100.txt"),  paste0(sam, "_viruses_names.txt")))
  }else{

    system(paste("STAR --runThreadN 8 --runMode alignReads --genomeDir",
                 "HumanGenome/ --sjdbGTFfile",
                 "HumanGenome/Genome.gtf --readFilesIn",
                 paste0("./fastq/", sam, "_1.fastq"),
                 "--outFilterMismatchNmax 2 --outFileNamePrefix", sam,
                 "--outReadsUnmapped Fastx --outSAMtype BAM SortedByCoordinate"))
    system(paste("rm -rf", paste0(sam, "\Aligned.sortedByCoord.out.bam"), paste0(sam, "\_STARgenome")))

    system(paste("bwa mem -t 8", "ViralGenome/viral_genome.fa",
                 paste0(sam, "Unmapped.out.mate1"), ">", paste0(sam, "_aln.sam")))

    system(paste("samtools view -Sb -h", paste0(sam, "_aln.sam"), ">", paste0(sam, "_aln.bam")))

    system(paste("samtools view", paste0(sam, "_aln.bam"), "| cut -f3 | sort | uniq -c | awk '{if ($1>=10000) print $0}' >", paste0(sam, "_viruses_count.txt")))

    system(paste("samtools sort", paste0(sam, "_aln.bam"), ">", paste0(sam, "_aln_sorted.bam")))

    system(paste("samtools depth", paste0(sam, "_aln_sorted.bam"), "| awk '{if ($3>=5) print $0}' | awk '{ if ($2!=(ploc+1)) {if (ploc!=0){printf(\"%s %d-%d\\n\",$1,s,ploc);}s=$2} ploc=$2; }' >", paste0(sam, "_continuous_region.txt")))

    system(paste("cut -d ' ' -f2", paste0(sam, "_continuous_region.txt"), "| awk -F '[-]' '{print$1\"\\t\"$2\"\\t\"$2-$1}' >", paste0(sam, "_ids.txt")))

    system(paste("paste", paste0(sam,"_continuous_region.txt"), paste0(sam,"_ids.txt"),"|awk '{if($5 >=100)print$1\"\\t\"$3\"\\t\"$4\"\\t\"$5}' >", paste0(sam,"_continuous_region_100.txt")))

    system(paste("awk '{print $2}'", paste0(sam, "_viruses_count.txt"), ">",   paste0(sam, "_viruses_names.txt")))

    system(paste("fgrep -wf", paste0(sam, "_viruses_names.txt"), paste0(sam,"_continuous_region_100.txt"), ">", paste0(sam,"_final_virus_list.txt")))

    system(paste("rm", paste0(sam,"_ids.txt"), paste0(sam,"_continuous_region_100.txt"), paste0(sam,"_viruses_names.txt")))
  }
}




