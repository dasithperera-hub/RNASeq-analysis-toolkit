# RNASeq-analysis-toolkit beta version
This repository contains (hopefully) an easy to use tool to analyze RNASeq data. 

## Dependencies
bowtie2 http://bowtie-bio.sourceforge.net/bowtie2/index.shtml

HTSeq https://htseq.readthedocs.io/en/master/


R packages

DESeq2 https://bioconductor.org/packages/release/bioc/html/DESeq2.html

plyr

dplyr

splitstackshape

ggplot2

RColorBrewer

EnhancedVolcano


The tool assumes all pre-processing steps have been carried out. This includes:

1) Quality control

Sequence quality can be assessed using tools such as fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Thereafter reads can be trimmed using tools such as trimmomatic http://www.usadellab.org/cms/?page=trimmomatic or cutadapt https://cutadapt.readthedocs.io/en/stable/

2) Assessing strandedness 

A script I haved used before can be found here: http://rseqc.sourceforge.net/#infer-experiment-py

3) Genome indexes have been built using a tool such as bowtie build http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

4) Genome annotation files (GFF) are available. 

These can be obtained from NCBI or rast can be used to generate a GFF file. https://rast.nmpdr.org/ 

5) Kegg ko numbers for the species of interest has been obtained.

This can be done using tools like BlastKOALA or GhostKOALA https://www.kegg.jp/blastkoala/

A parsed, up to date version of the kegg ko database would be ideal https://www.genome.jp/kegg-bin/get_htext?ko00001. A parsed version is uploaded here along with the scripts. 




## How to use- important information current version

1) The fastq files require suffixes in the following format:  1_R1.fastq 1_R2.fastq, where the first digit is the rep number and R1 or R2 is the read number. 

2) The file containing ko numbers for the species needs to have the suffix:  _ko.txt

3) The kegg database file needs to have the suffix: _kegg.csv

4) The current version assumes that the fastq files, gff file, kegg database file and ko file are in the working directory.

5) The Shell script and R script need to be in working directory.

6) The current version assumes the index and gff files have the same name. for e.g. if index was ecoli, gff file should be named ecoli.gff 

7) If HTSeq generated count files are already available the analysis can be carried out using standaloneR.sh. The count files would need the prefix "condition1_" or "condition2_" and the suffix "_counts#.txt", where # is the number of the replicate. For eg. the third replicate for the second condition called ecan should have a count file named "condition2_ecan_counts3.txt".       


## Usage
 
sbatch paired_end_duplicate_analysis.sh -d (directory/path to fastq files) -1 (prefix of condition 1 fastq files) -2 (prefix of condition 1 fastq files) -i (index file) -p (path to index files) -s (strandedness) -f (feature). 


## Usage example

Following files are in working directory:

GFF file: ecoli.gff

Kegg database table : May2020_kegg.csv

Ecoli ko table: ec_ko.txt

deseq_tablesandkegg.R

paired_end_duplicate_analysis.sh


Not in directory:

bowtie index: ecoli

path to bowtie index is /data/mramseylab/genome/ecoli/

Condition 1 fastq files: ec1_R1.fastq ec1_R2.fastq ec2_R1.fastq ec2_R2.fastq

Condition 2 fastq files: ecan1_R1.fastq ecan1_R2.fastq ecan2_R1.fastq ecan2_R2.fastq

path to fastq files: /data/mramseylab/raw_reads/anaerobic/



The intention is to check for DEG coding sequences using the reverse strandedness parameter.

Run the following script:

sbatch paired_end_duplicate_analysis.sh -d /data/mramseylab/raw_reads/anaerobic/ -1 ec -2 ecan -i ecoli -p /data/mramseylab/genome/ecoli/ -s reverse -f CDS
