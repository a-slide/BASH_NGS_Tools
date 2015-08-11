#! /bin/bash

if [[ "$#" != 7 ]]
then
    echo -e "USAGE: RNAseq_U7_pipeline.sh <STAR_index_dir> <R1.fa.gz> <R2.fa.gz> <gff_file> <out_prefix> <min_mapq_score> <n_threads>"
    exit 1
fi

# OPTIONS
# $1 = STAR_index_dir
# $2 = R1.fa.gz
# $3 = R2.fa.gz
# $4 = gff_file
# $5 = out_prefix
# $6 = min_mapq_score
# $7 = n_threads

# DEPENDENCIES
# samtools 1.2+
# STAR 2.4.2.a+
# samstat 1.5+ 
# bedtools 2.17.0+
# HTSeq 0.6.1p1+

echo -e "\n########################## PROCESSING SAMPLE: $5 ##########################\n"

# Add a dash to the file prefix 
bn="$5"_

# Align with STAR
echo -e "Mapping with STAR and sort bam\n"
STAR --runThreadN "$7" --readFilesCommand "gunzip -c" --genomeDir "$1" --readFilesIn "$2" "$3" \
--outFileNamePrefix "$bn" --quantMode GeneCounts --outSAMtype BAM SortedByCoordinate \
--outFilterType BySJout

echo -e "\nGenerate a report\n"
samstat "$bn"Aligned.sortedByCoord.out.bam

## Filter out unmapped and secondary reads 
echo -e "Filter_bam\n"
samtools view "$bn"Aligned.sortedByCoord.out.bam -hb -F 260 -q "$6" -@ "$7" | samtools sort - -@ "$7" -o a > "$bn".bam

## Index bam and generates bedgraph
echo -e "\nCreating a bam index\n"
samtools index "$bn".bam &

echo -e "\nCreating bedgraphs\n"
genomeCoverageBed -ibam "$bn".bam -bg -split -strand + -trackline -trackopts "name="$bn"Positive_Strand color=250,0,0" | \
awk '$4 >= 5' > "$bn"positive_strand.bedgraph &
genomeCoverageBed -ibam "$bn".bam -bg -split -strand - -trackline -trackopts "name="$bn"Negative_Strand color=250,0,0" | \
awk '$4 >= 5' > "$bn"negative_strand.bedgraph

# Create HTSeq-Count file
echo -e "\nCounting reads overlaping features with HTSeq count\n"
htseq-count -f bam -s yes -m union -t exon "$bn".bam "$4" > "$bn".counts

echo -e "\n### DONE ###\n"

exit 0
