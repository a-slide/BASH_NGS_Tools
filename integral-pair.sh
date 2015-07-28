#! /bin/bash

	###################################################################################
	#
	#
	#
	###################################################################################
	
	if [[ $# == 6 ]]; then
	
		#Tri préalable sous des fastq 
		
		echo -e "\n\tTrimming bad quality sequences...\n" 
		sickle pe -f "$4" -r "$5" -t sanger -o "$6"_trim_R1.fastq -p "$6"_trim_R2.fastq -s "$6"_trim_single.fastq -q 30 -l 30
		
		# Mapping against reference 1 index and extracting of fully unmapped paired reads
		
		mkdir -p step1
		cd step1
		echo -e "\n\tProcessing raw sequences against reference 1 index...\n" 
		extract-unmap.sh ../$1 ../"$6"_trim_R1.fastq ../"$6"_trim_R2.fastq "$6"_1
		cd ..
		
		# Second step of mapping against reference 2 index and extracting of fully unmapped paired reads
		
		mkdir -p step2
		cd step2
		echo -e "\n\tProcessing unmapped fastq from the first step of mapping against reference 2 index...\n" 
		extract-unmap.sh ../"$2" ../step1/"$6"_1_R1.fastq ../step1/"$6"_1_R2.fastq "$6"_2
		cd ..
		
		# Last step of mapping against a mixed ref1 ref2 index
		
		mkdir -p step3
		cd step3
		
		echo -e "\n\tProcessing unmapped fastq from the second step of mapping against reference 1+2 index...\n" 
		
		echo -e "\n\tLocal mapping with Bowtie2...\n" 
		if [[ -f "$6"_3.sam ]]
			then echo -e "\n\tOutput file already existing. Skiping step\n" 
		else
			bowtie2 --very-sensitive-local -p 12 --reorder -x ../"$3" -1 ../step2/"$6"_2_R1.fastq -2 ../step2/"$6"_2_R2.fastq -S "$6"_3.sam
		fi
		
		echo -e "\n\tConverting SAM output in sorted BAM file with Samtools...\n"
		if [[ -f "$6"_3.bam ]]
			then echo -e "\n\tOutput file already existing. Skiping step\n" 
		else
			samtools view -hbS "$6"_3.sam | samtools sort - "$6"_3
		fi
		
		echo -e "\n\tRemoving reads unmapped and with  MAPQ < 5...\n"
		if [[ -f "$6"_3_mapped.bam ]]
			then echo -e "\n\tOutput file already existing. Skiping step\n" 
		else
			samtools view -h -F 4 -q 5 -b "$6"_3.bam > "$6"_3_mapped.bam
		fi
		
		echo -e "\n\tRemoving PCR duplicates...\n"
		if [[ -f "$6"_3_mapped_rmdup.bam ]]
			then echo -e "\n\tOutput file already existing. Skiping step\n" 
		else
			samtools rmdup -s "$6"_3_mapped.bam "$6"_3_mapped_rmdup.bam
		fi
		
		echo -e "\n\tStarting convertion in SAM file with samtools..."
		if [[ -f "$6"_3_mapped_rmdup.sam ]]
			then echo -e "\n\tOutput file already existing. Skiping step\n" 
		else
			samtools view -h "$6"_3_mapped_rmdup.bam > "$6"_3_mapped_rmdup.sam
		fi
		
	# Creation d'un bedgraph pour visualisation sur UCSC
		
		echo -e "\n\tStarting generating BedGraph..." 
		if [[ -f "$6".bedGraph ]]
			then echo -e "\n\tOutput file already existing. Skiping step\n" 
		else
			genomeCoverageBed -ibam "$6"_3_mapped_rmdup.bam -bg -trackline -trackopts "name="$6" color=250,0,0" > "$6".bedGraph
		fi
		
		cd ..
		
	# Creation d'un rapport contenant les flagstats
		echo -e "\n\tGenerating a flagstat report"
		
		if [[ -f "$6"_report.txt ]]
			then echo -e "\n\tOutput file already existing. Skiping step\n" 
		else
			echo -e "\n\tReference 1 path = "$1"" > "$6"_report.txt
			echo -e "\n\tRaw mapping\n" >> "$6"_report.txt
			samtools flagstat step1/"$6"_1.bam  >> "$6"_report.txt
			echo -e "\n\tUnmapped sequences\n" >> "$6"_report.txt
			samtools flagstat step1/"$6"_1_unmapped.bam  >> "$6"_report.txt
			
			echo -e "\n\tReference 2 path = "$2"" >> "$6"_report.txt
			echo -e "\n\tRaw mapping\n" >> "$6"_report.txt
			samtools flagstat step2/"$6"_2.bam  >> "$6"_report.txt
			echo -e "\n\tUnmapped sequences\n" >> "$6"_report.txt
			samtools flagstat step2/"$6"_2_unmapped.bam  >> "$6"_report.txt
			
			echo -e "\n\tReference 1+2 path = "$3"" >> "$6"_report.txt
			echo -e "\n\tRaw mapping\n" >> "$6"_report.txt
			samtools flagstat step3/"$6"_3.bam >> "$6"_report.txt
			echo -e "\n\tFinal result with removal of PCR duplicates\n" >> "$6"_report.txt
			samtools flagstat step3/"$6"_3_mapped_rmdup.bam >> "$6"_report.txt
		
		# Comptage du nombre de hits dans chaque sequence de reférence 
		#	echo -e "\n\tAdding a summary of hits to the report"
		#	echo -e "\n\tSummary of hits\n" >> "$6"_report.txt
		#		for hit in `cat "$2".list`
		#			do
		#			n=$(egrep -c ""$hit"\s" "$5"_host_mapped_rmdup.sam)
		#			n=$((--n))
		#			echo "$hit : "$n"" >> "$5"_report.txt
		#			#egrep -c ""$hit"\s" "$5"_host_mapped_rmdup.sam >> "$5"_report.txt
		#		done
		
			echo -e "\n\tPrinting the report\n"
			cat "$6"_report.txt
		fi
	
	# rangement dossier
	#	gzip step*/*.fastq
	
		exit 0
	
	else 
		echo -e "\t6 parameters are required \n"
		echo -e "\t1 = Reference 1 index-file"
		echo -e "\t2 = Reference 2 index-file"
		echo -e "\t3 = Reference 1+2 index-file"
		echo -e "\t4 = R1.fastq(.gz)"
		echo -e "\t5 = R2.fastq(.gz)"
		echo -e "\t6 = Output.name\n"
		exit 1
	fi
