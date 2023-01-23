# make a new directory called genome 
mkdir HumanGenome
cd HumanGenome

# download the "Genome sequence, primary assembly (GRCh38)" fasta file
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/GRCh38.primary_assembly.genome.fa.gz
# filter it as described in the note below
# download the annotations that correspond to it 
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_29/gencode.v29.annotation.gtf.gz

gunzip GRCh38.primary_assembly.genome.fa.gz
mv GRCh38.primary_assembly.genome.fa Genome.fa
# wait
gunzip gencode.v29.annotation.gtf.gz
mv gencode.v29.annotation.gtf Genome.gtf


GENOMEDIR="HumanGenome"
mdkir -p $GENOMEDIR/STAR
STAR --runThreadN 23 --runMode genomeGenerate --genomeDir $GENOMEDIR/STAR --genomeFastaFiles $GENOMEDIR/GRCh38.primary_assembly.genome.fa --sjdbGTFfile $GENOMEDIR/gencode.v29.annotation.gtf


cd ../
mkdir ViralGenome
wget https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.1.genomic.fna.gz
gunzip viral.1.1.genomic.fna.gz
mv viral.1.1.genomic.fna viral_genome.fa
bwa index index -a bwtsw viral_genome.fa viral_genome.fa