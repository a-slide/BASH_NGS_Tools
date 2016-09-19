#! /bin/bash

if [[ "$#" != 6 ]]
then
    echo -e "USAGE: RNAseq_U7_pipeline.sh <STAR_index_dir> <R1.fa.gz> <R2.fa.gz> <gff_file> <out_prefix> <min_mapq_score> <n_threads>"
    exit 1
fi

# OPTIONS
# $1 = STAR_index_dir
# $2 = R1.fa.gz
# $3 = R2.fa.gz
# $4 = n_threads
# $5 = out_dir
# $6 = min_mapq_score

# DEPENDENCIES
# samtools 1.2+
# STAR 2.4.2.a+
# samstat 1.5+ 
# bedtools 2.17.0+

echo -e "\n########################## PROCESSING SAMPLE: $5 ##########################\n"

mkdir $5

# Align with STAR
echo -e "Mapping with STAR and sort bam\n"
STAR \
	--genomeDir "$1"\
	--readFilesIn "$2" "$3" \
	--outFileNamePrefix ./"$5"/ \
	--runThreadN "$4"\
	--readFilesCommand "gunzip -c"\
	--quantMode GeneCounts\
	--outSAMtype BAM SortedByCoordinate \
	--outFilterMultimapNmax 20\
	--alignSJoverhangMin 8\
	--alignSJDBoverhangMin 1\
	--outFilterMismatchNmax 999\
	--alignIntronMin 20\
	--alignIntronMax 1000000\
	--alignMatesGapMax 1000000\
	--outFilterIntronMotifs RemoveNoncanonical\
	--outFilterScoreMin "$6"

# Rename the output bam file
mv ./"$5"/Aligned.sortedByCoord.out.bam ./"$5"/out.bam
 
## Generate a report, index bam and generates bedgraphs
echo -e "\nGenerate a report\n"
samstat ./"$5"/out.bam &
echo -e "\nCreating a bam index\n"
samtools index ./"$5"/out.bam &
echo -e "\nCreating bedgraphs\n"
genomeCoverageBed -ibam ./"$5"/out.bam -bg -split -strand + -trackline -trackopts "name="$5"_Positive_Strand color=250,0,0" | \
awk '$4 >= 5' > ./"$5"/positive_strand.bedgraph &
genomeCoverageBed -ibam ./"$5"/out.bam -bg -split -strand - -trackline -trackopts "name="$5"_Negative_Strand color=250,0,0" | \
awk '$4 >= 5' > ./"$5"/negative_strand.bedgraph

echo -e "\n### DONE ###\n"

exit 0
