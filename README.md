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


## Usage
 
sbatch paired_end_duplicate_analysis.sh -1 (prefix of condition 1 fastq files) -2 (prefix of condition 1 fastq files) -i (index file) -p (path to index file) -s (strandedness) -f (feature). 






 


