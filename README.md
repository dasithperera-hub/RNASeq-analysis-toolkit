# RNASeq-analysis-toolkit beta version
This repository contains (hopefully) an easy to use tool to analyze RNASeq data. The tool assumes all pre-processing steps have been carried out. This includes:

1) Quality control
Sequence quality can be assessed using tools such as fastqc https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

Thereafter reads can be trimmed using tools such as trimmomatic http://www.usadellab.org/cms/?page=trimmomatic or cutadapt https://cutadapt.readthedocs.io/en/stable/

2) Assessing strandedness 
A script I haved used before can be found here: http://rseqc.sourceforge.net/#infer-experiment-py

3) Genome indexes have been built using a tool such as bowtie build http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml

4) Genome annotation files (GFF) are available. These can be obtained from NCBI or rast can be used to generate a GFF file. https://rast.nmpdr.org/ 

5) Kegg ko numbers for the species of interest has been obtained
This can be done using tools like BlastKOALA or GhostKOALA https://www.kegg.jp/blastkoala/

A parsed, up to date version of the kegg ko database would be ideal https://www.genome.jp/kegg-bin/get_htext?ko00001. A parsed version is uploaded here along with the scripts. 


## How to use- important information beta version




 


