# ViralCellDetector
This R script can be used to detect virus infection in cell lines. It takes the input as raw RNA-seq reads as input and provide the list of putative viruses as output. The final virus could be visualized by the igv genome viewer. 

## Installation
```
git clone https://github.com/Bin-Chen-Lab/ViralCellDetector

## Then go to the ViralCellDetector directory:
cd ViralCellDetector
```
# Usage
First user need to run the Genome_index.sh bash code
```r
bash Genome_index.sh
## This will download the Human genome and Viral genomes required for the script.
````

To use the script, you will need to have R and the necessary packages installed. Once you have the dependencies set up, you can run the script by calling following code
``` r
## basic example code
Rscript ViralCellDetector.R sample_input.txt
## The sample_inpit.txt file is the file containing the name of fastq input files. If a fastq input file is "input_file_1.fq and input_file_2.fq", then "sample_input.txt" file should have "input_file" name in first row.
```

# Input

The script takes as input as RNA-seq raw sequencing fastq file, Users need to create a fastq directory in the ViralCellDector directory and put all their fastq files in the fastq directory. Then user need to creat a sample_input.txt file and enter all fastq files name in it. If fastq files are paired end then their only name should be provided in sample_input.txt as mentioned in example code. The ViralCellDetector prefer to accept the RNA-seq data  generated by the robominus approach, however user can also use the RNA-seq raw reads generated by poly-A enrichment approach.

# Output

The script will output a list of names of putative viruses, which are supposed to be present in the samples. 

# Note

Please make sure that the sequencing reads are with good quality.

# Dependencies

The script requires the following R packages and tools to be installed:
samtools (https://github.com/samtools/samtools/)

STAR (https://github.com/alexdobin/STAR)

BWA (https://bio-bwa.sourceforge.net/)

dplyr (https://dplyr.tidyverse.org/)

plyr (https://cran.r-project.org/web/packages/plyr/index.html)



Make sure these packages are installed before running the script.

# Example

To run the example files, we have provided a example code:

```r
Run_test_ViralCellDetector.sh
```
# Contributors

ViralCellDetector is developed by BinChen lab. Any questions or feedback can be
addressed to Rama Shankar, PhD, <ramashan@msu.edu> or Bin Chen, PhD, PI,
<Chenbi12@msu.edu>


