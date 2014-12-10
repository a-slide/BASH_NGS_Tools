#! /bin/bash


if [[ $# < 3 ]]; then
    echo -e "\nUsage of HTS_mapping_AAV_and_back\n
    Arg 1 = index
    Arg 2 = Sample1_R1.fastq
    Arg 3 = Suffix of R1 to remove from outpit name
    Arg 4 = Sample1_R2.fastq
    [Arg n-1 = SampleN_R1.fastq
    Arg n = SampleN_R2.fastq]"
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
    echo -e "Mapping with bwa mem and sort bam\n"
    bwa-0.7.10 mem -M -t 12 "$INDEX" "$1" "$2" | samtools view -F 4 -hbS - | samtools sort - "$OUTPUT"
    echo -e "\nProcessing with samtools index\n"
    samtools index "$OUTPUT".bam
    echo -e "\nCreating bed file\n"
    genomeCoverageBed -ibam "$OUTPUT".bam -bga -trackline -trackopts "name="$OUTPUT" color=250,0,0" > "$OUTPUT".bedgraph
    # Shift the arguments to process next argument 
    shift
    shift
done

exit 0
