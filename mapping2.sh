#! /bin/zsh

	if [[ $# == 5 ]]; then
		echo "Starting mapping with bowtie 2 \n " 
		
# Mapping avec bowtie2 en utilisant différentes options
#  -local = alignement local avec soft cliping
# -p 4 = 4 thread en parallèle pour proc multithread option --reorder pour remettre les sequences dans l'ordre initial
# --no-hd pour supprimer le header du fichier SAM
 
 
		bowtie2 --local -p 4 --reorder -x "$1" -1 "$2"  -2 "$3" -S "$4".sam
		
	else 
		echo " \n5 parameters required. Use the following synthax \n"
		echo "mapping.sh [index-file] [R1.gz] [R2.gz] [output.name] [text file with list of fasta headers]\n "
		exit 1	
	fi
	
	echo " \n \nStarting convertion in sorted BAM file with sam tools \n "
	
	samtools view -bS "$4".sam | samtools sort - "$4"
	
	echo " \n \nSimple stats using SAMTools flagstat \n "
	samtools flagstat "$4".bam
	samtools flagstat "$4".bam > "$4"_flagstat.txt
	
	echo " \n \nRemoving unmapped reads \n "
	
	samtools view -h -F 4 -b "$4".bam > "$4"_mapped.bam
	
	echo " \n \nSimple stats using SAMTools flagstat \n "
	samtools flagstat "$4"_mapped.bam
	samtools flagstat "$4"_mapped.bam > "$4"_mapped_flagstat.txt
	
#Creation d'un fichier qui recapitule le nombre de hits dans chaque sequence de reférence 
#trouver un moyen de virer le header du sam = 95 premières lignes
	
	echo "Summary of hits \n" > "$4"_hit.txt

# Pour compter les hits dans les sequences de ref

	for hit in `cat $5`
	do
		echo "$hit : " >> "$4"_hit.txt
		cat "$4".sam | egrep -c "$hit" >> "$4"_hit.txt
	done
	
	less "$4"_hit.txt

	exit 0
