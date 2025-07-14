# ViralCellDetector

**ViralCellDetector** is an R-based pipeline designed to detect viral infections in RNA-seq samples. It accepts raw RNA-seq FASTQ files as input and outputs a list of putative viruses identified in the sample. The final results can be visualized using genome browsers such as IGV (Integrative Genomics Viewer).

---

## 📦 Installation

```bash
git clone https://github.com/Bin-Chen-Lab/ViralCellDetector
cd ViralCellDetector
```

---

## 🚀 Usage

### Step 1: Prepare the Reference Genomes

Before running the detection script, you need to download the host genome and annotation files by executing the provided shell script.

1. Edit the `Genome_file.txt` to include the appropriate FTP links for your species' genome and annotation files. For example, for the human genome:

```
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz  
ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz
```

2. Then run:

```bash
bash Genome_index.sh Genome_file.txt
```

This will:
- Download the host genome and annotation files
- Download the viral genomes from NCBI
- Build STAR and BWA index files for alignment

---

### Step 2: Prepare Input Files

1. Create a folder named `fastq/` inside the `ViralCellDetector` directory and move your FASTQ files into this folder.

2. Create a file named `sample_input.txt` containing the **base names** (without `_1.fastq` or `_2.fastq` suffixes) of your paired-end FASTQ files.  
   Example entry for a sample:
   ```
   input_file
   ```

If the paired-end files are `input_file_1.fastq` and `input_file_2.fastq`, only write `input_file` in the `sample_input.txt`.

---

### Step 3: Run the Pipeline

After the setup, execute the main R script:

```r
Rscript ViralCellDetector.R sample_input.txt
```

---

## 🧬 Input

- Raw RNA-seq FASTQ files (preferably rRNA-depleted, though poly-A data can also be used)
- A text file (`sample_input.txt`) listing sample base names
- The `fastq/` directory should contain all input FASTQ files

---

## 📤 Output

- A summary file for each sample listing:
  - Putative virus names
  - Genome size
  - Number of mapped reads
  - Covered genomic regions

These results can be loaded into IGV or other genome browsers for visualization.

---

## ⚠️ Notes

- Ensure high-quality sequencing reads for optimal detection.
- Poly-A enriched data may result in reduced viral diversity compared to rRNA-depleted data.

---

## 🧩 Dependencies

Make sure the following tools and R packages are installed before running the script:

### 🔧 Command-line Tools
- [`samtools`](https://github.com/samtools/samtools)
- [`STAR`](https://github.com/alexdobin/STAR)
- [`BWA`](https://bio-bwa.sourceforge.net/)

### 📦 R Packages

Install these via CRAN if not already available:

```r
install.packages(c("dplyr", "data.table", "tidyr"))
```

---

## 🧪 Example

To run the example test case:

```bash
bash Run_test_ViralCellDetector.sh
```

---

## 👥 Contributors

**ViralCellDetector** is developed and maintained by the [Bin Chen Lab](https://binchenlab.org/).

For questions or suggestions, please contact:

- Rama Shankar, PhD – <ramashan@msu.edu>  
- Bin Chen, PhD (PI) – <Chenbi12@msu.edu>
