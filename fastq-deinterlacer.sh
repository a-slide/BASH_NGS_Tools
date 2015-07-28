#! /bin/bash

if [[ $# == 2 ]]; then

	sed -n '1~8{N;N;N;p}' "$1" | gzip > "$2"_R1.fastq.gz
	sed -n '5~8{N;N;N;p}' "$1" | gzip > "$2"_R2.fastq.gz
	exit 0
	
else
	echo -e "Wrong number of arguments.\n"
	echo -e "fastq-desinterlacer [sorted fastq file with paired ends to deinterlace] [output name]"
	exit 1
fi
