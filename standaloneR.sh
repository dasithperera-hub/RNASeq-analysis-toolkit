#!/bin/bash
#SBATCH -t 0:30:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE


module load R-bundle-Bioconductor/3.10-foss-2019b

mkdir figures
mkdir kegg_analysis
mkdir rnaseq_analysis

Rscript deseq_tablesandkegg.R
