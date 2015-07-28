#! /bin/bash

if [[ "$#" != 4 ]]
then
    echo -e "USAGE: bwa_pipeline_single.sh n_thread index_path fasta_path out_suffix"
    exit 1
fi

echo -e "Mapping with bwa mem and sort bam\n"
bwa mem -M -t "$1" "$2" "$3" | samtools view -F 4 -hbS - | samtools sort - "$4"

echo -e "\nProcessing with samtools index\n"
samtools index "$4".bam

echo -e "\nCreating bedgraph\n"
genomeCoverageBed -ibam "$4".bam -bga -trackline -trackopts "name="$4" color=250,0,0" > "$4".bedgraph

exit 0
