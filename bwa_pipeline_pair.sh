#! /bin/bash

if [[ "$#" == 0 ]]
then
    echo -e "\nUsage of HTS_mapping_AAV_and_back <INDEX> <Sufix_to_remove> <S1-R1.fastq> <S1-R2.fastq> [<S2-R1.fastq> <S2-R2.fastq>  ...]"
    exit 1
fi

# Save arguments 1 and 2 and shift them
INDEX=$1
SUFFIX=$2
shift
shift

# Define the range of argument to iterate through
START=1
END=$#

# Iterate over arguments by steps of 2
for ((i=START ; i<=END ; i+=2))
do
    
    # Generate an output name by extracting the basename and the suffix from first fastq file name
    b=`basename $1`
    OUTPUT=${b%$SUFFIX}
    
    echo -e "PROCESSING $OUTPUT : ALIGNING `basename $1` AND `basename $2` AGAINST `basename $INDEX`"
    # Main pipeline
    echo -e "Mapping with bwa mem and sort bam\n"
    bwa-0.7.10 mem -M -t 12 "$INDEX" "$1" "$2" | samtools view -F 4 -hbS - | samtools sort - "$OUTPUT"
    echo -e "\nProcessing with samtools index\n"
    samtools index "$OUTPUT".bam
    echo -e "\nCreating bed file\n"
    genomeCoverageBed -ibam "$OUTPUT".bam -bga -trackline -trackopts "name="$OUTPUT" color=250,0,0" > "$OUTPUT".bedgraph
    genomeCoverageBed -ibam "$OUTPUT".bam -d > "$OUTPUT".bed
    echo -e "$OUTPUT PROCESSED\n"
    
    # Shift the arguments to process next fastq files
    shift
    shift
done

echo -e "#### DONE ####"

exit 0
