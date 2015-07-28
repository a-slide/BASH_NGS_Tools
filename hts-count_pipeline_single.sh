#! /bin/bash

if [[ "$#" < 3 ]]
then
    echo -e "USAGE: HTS-count_pipeline_single.sh <gff_file> <bam_file1> [<bam_file2> <bam_file3>...]"
    exit 1
fi

GFF=$1
shift

for i in $@
do
    echo -e "Counting reads per gene in $i"
    bsub -q research-rh6 -M 20000 -R 'rusage[mem=20000] select[mem>20000]' -n 1 -e log.stderr -o log.stdout "htseq-count -f bam $i $GFF -s no -m intersection-nonempty > ./$i.counts"

done

exit 0
