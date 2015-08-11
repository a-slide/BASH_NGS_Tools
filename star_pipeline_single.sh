#! /bin/bash

if [[ "$#" != 4 ]]
then
    echo -e "USAGE: star_pipeline_single.sh n_thread index_path fasta_path out_suffix"
    exit 1
fi

echo -e "Mapping with star and sort bam\n"

STAR --runThreadN "$1" --genomeDir "$2" --readFilesIn "$3" --outFileNamePrefix "$4"
samtools view -F 4 -hbS "$4"Aligned.out.sam | samtools sort - "$4"
rm "$4"Aligned.out.sam

echo -e "\nProcessing with samtools index\n"
samtools index "$4".bam

echo -e "\nCreating bedgraph\n"
genomeCoverageBed -ibam "$4".bam -bga -trackline -trackopts "name="$4" color=250,0,0" > "$4".bedgraph

exit 0
