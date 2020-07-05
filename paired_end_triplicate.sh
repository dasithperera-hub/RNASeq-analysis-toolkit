#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --nodes=1 --ntasks-per-node=1
#SBATCH --export=NONE

module load Bowtie2/2.3.5.1-GCC-8.3.0
module load HTSeq/0.11.2-foss-2019b-Python-3.7.4
module load R-bundle-Bioconductor/3.10-foss-2019b


while getopts "d:1:2:i:p:s:f:" opt; do
  case $opt in
    d) DIRFASTQ=$OPTARG     ;;  
    1) CONDITION1=$OPTARG      ;;
    2) CONDITION2=$OPTARG   ;;
    i) INDEX=$OPTARG	;;
    p) DIRPATH=$OPTARG	;;
    s) STRANDED=$OPTARG ;;
    f) FEATURE=$OPTARG ;;
  esac
done

if [ -z "$DIRFASTQ" ]; then
    echo "directory/path to fastqs were NOT given, exit."
    exit 1;
fi

if [ -z "$CONDITION1" ]; then
    echo "condition1 was NOT given, exit."
    exit 1;
fi

if [ -z "$CONDITION2" ]; then
    echo "condition2 was NOT given, exit."
    exit 1;
fi

if [ -z "$INDEX" ]; then
    echo "index was NOT given, exit."
    exit 1;
fi

if [ -z "$DIRPATH" ]; then
    echo "path to index was NOT given, exit."
    exit 1;
fi

if [ -z "$STRANDED" ]; then
    echo "stranded parameter was NOT given, exit."
    exit 1;
fi


bowtie2 -x $DIRPATH/$INDEX -1 $DIRFASTQ"$CONDITION1"1_R1.fastq -2 $DIRFASTQ"$CONDITION1"1_R2.fastq -S "$CONDITION1"_SAMPLES1.sam
bowtie2 -x $DIRPATH/$INDEX -1 $DIRFASTQ"$CONDITION1"2_R1.fastq -2 $DIRFASTQ"$CONDITION1"2_R2.fastq -S "$CONDITION1"_SAMPLES2.sam
bowtie2 -x $DIRPATH/$INDEX -1 $DIRFASTQ"$CONDITION1"3_R1.fastq -2 $DIRFASTQ"$CONDITION1"3_R2.fastq -S "$CONDITION1"_SAMPLES3.sam

bowtie2 -x $DIRPATH/$INDEX -1 $DIRFASTQ"$CONDITION2"1_R1.fastq -2 $DIRFASTQ"$CONDITION2"1_R2.fastq -S "$CONDITION2"_SAMPLES1.sam
bowtie2 -x $DIRPATH/$INDEX -1 $DIRFASTQ"$CONDITION2"2_R1.fastq -2 $DIRFASTQ"$CONDITION2"2_R2.fastq -S "$CONDITION2"_SAMPLES2.sam
bowtie2 -x $DIRPATH/$INDEX -1 $DIRFASTQ"$CONDITION2"3_R1.fastq -2 $DIRFASTQ"$CONDITION2"3_R2.fastq -S "$CONDITION2"_SAMPLES3.sam

htseq-count -s $STRANDED -t $FEATURE -a 1 -m intersection-nonempty -i ID "$CONDITION1"_SAMPLES1.sam "$INDEX".gff > condition1_"$CONDITION1"_counts1.txt
htseq-count -s $STRANDED -t $FEATURE -a 1 -m intersection-nonempty -i ID "$CONDITION1"_SAMPLES2.sam "$INDEX".gff > condition1_"$CONDITION1"_counts2.txt
htseq-count -s $STRANDED -t $FEATURE -a 1 -m intersection-nonempty -i ID "$CONDITION1"_SAMPLES3.sam "$INDEX".gff > condition1_"$CONDITION1"_counts3.txt

htseq-count -s $STRANDED -t $FEATURE -a 1 -m intersection-nonempty -i ID "$CONDITION2"_SAMPLES1.sam "$INDEX".gff > condition2_"$CONDITION2"_counts1.txt
htseq-count -s $STRANDED -t $FEATURE -a 1 -m intersection-nonempty -i ID "$CONDITION2"_SAMPLES2.sam "$INDEX".gff > condition2_"$CONDITION2"_counts2.txt
htseq-count -s $STRANDED -t $FEATURE -a 1 -m intersection-nonempty -i ID "$CONDITION2"_SAMPLES3.sam "$INDEX".gff > condition2_"$CONDITION2"_counts3.txt

mkdir alignment
mv *.sam ./alignment

mkdir figures
mkdir kegg_analysis
mkdir rnaseq_analysis

Rscript deseq_tablesandkegg.R
