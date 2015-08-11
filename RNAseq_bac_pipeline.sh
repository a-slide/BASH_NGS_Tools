#! /bin/bash

if [[ "$#" != 5 ]]
then
    echo -e "USAGE: RNAseq_bac_pipeline.sh <bwa_index_path> <sequence.fa.gz> <gff_file> <min_mapq_score> <n_threads>"
    exit 1
fi

# OPTIONS
# $1 = STAR_index_dir
# $2 = R1.fa.gz
# $3 = gff_file
# $4 = min_mapq_score
# $5 = n_threads

# Dependencies
# samtools 1.2+
# bwa 0.1.12+
# samstat 1.5+ 
# bedtools 2.17.0+
# HTSeq 0.6.1p1+

# Extract basename
bn=`basename ${2%.fastq.gz}`

echo -e "\n########################## PROCESSING SAMPLE: $5 ##########################\n"

# Align with bwa mem
echo -e "Mapping with bwa mem and sort bam\n"
bwa mem -M -t "$5" "$1" "$2" | samtools view - -hb -@ "$5" | samtools sort - -@ "$5" -o a > "$bn"_raw.bam
samstat "$bn"_raw.bam

## Filter out unmapped and secondary reads 
echo -e "Filter_bam\n"
samtools view "$bn"_raw.bam -hb -F 260 -q "$4" -@ "$5"| samtools sort - -@ "$5" -o a > "$bn".bam
samstat "$bn".bam 

## Index bam and generates bedgraph
echo -e "\nProcessing with samtools index\n"
samtools index "$bn".bam
echo -e "\nCreating bedgraph\n"
genomeCoverageBed -ibam "$bn".bam -bg -trackline -trackopts "name="$bn" color=250,0,0" | \
awk '$1 ~ "bedGraph" || $4 > 5' > "$bn".bedgraph

# Create HTSeq-Count file
echo -e "\nCounting reads overlaping features with HTSeq count\n"
htseq-count -f bam -s no -m union -t exon "$bn".bam "$3" > "$bn".counts

exit 0
