#! /bin/bash

	if [[ $# == 1 ]]; then
		
		echo "Generation of the list"
		egrep "^>" $1.fa | sed 's/>//g' > "$1".list
		exit 0
		
	else 
		echo "Wrong number of arguments"
		echo "filelist.sh [fasta path without extension]"
		exit 1	
	
	fi
