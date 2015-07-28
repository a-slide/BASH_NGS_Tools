#! /bin/bash

	###################################################################################
	#
	#
	#
	###################################################################################
	
	if [[ $# == 6 ]]; then
	
		echo -e "\n\tPre-processing of reads...\n" 
	
		#Fusion des read pair-ends 
		
		echo -e "\nMerging R1 and R2...\n" 
		if [[ -f "$6"_merged.fastq ]]
			then echo -e "Output file already existing. Skiping step\n" 
		else
			gzip -d -c "$4" > "$6"_merged.fastq
			gzip -d -c "$5" >> "$6"_merged.fastq
		fi
		
		#Quality filter préalable avec Sickle
		
		echo -e "\nRemoving bad quality sequences (q 25 -l 151)...\n" 
		if [[ -f "$6"_merged_fitered.fastq ]]
			then echo -e "Output file already existing. Skiping step\n" 
		else
			sickle se -t sanger -q 25 -l 151 -n -f "$6"_merged.fastq -o "$6"_merged_fitered.fastq
		fi
		
		# Mapping against reference 1 index and extracting of fully unmapped paired reads
		
		mkdir -p step1
		cd step1
		echo -e "\n\tProcessing raw sequences against reference 1 index...\n" 
		extract-unmap-single.sh ../$1 ../"$6"_merged_fitered.fastq "$6"_1
		cd ..
		
		# Second step of mapping against reference 2 index and extracting of fully unmapped paired reads
		
		mkdir -p step2
		cd step2
		echo -e "\n\tProcessing unmapped fastq from the first step of mapping against reference 2 index...\n" 
		extract-unmap-single.sh ../"$2" ../step1/"$6"_1_unmapped.fastq "$6"_2
		cd ..
		
		# Last step of mapping with BWA against virus index
		
		mkdir -p step3
		cd step3
		
		echo -e "\n\tProcessing unmapped fastq from the second step of mapping against reference 3 index...\n" 
		
		echo -e "\nLocal mapping with BWAsw...\n" 
		if [[ -f "$6"_3.sam ]]
			then echo -e "Output file already existing. Skiping step\n" 
		else
			bwa bwasw -t 12 ../"$3" ../step2/"$6"_2_unmapped.fastq >"$6"_3.sam
		fi
		
		echo -e "\nConverting SAM output in sorted BAM file with Samtools...\n"
		if [[ -f "$6"_3.bam ]]
			then echo -e "Output file already existing. Skiping step\n" 
		else
			samtools view -hbS "$6"_3.sam | samtools sort - "$6"_3
		fi
		
		echo -e "\nRemoving reads unmapped and with  MAPQ < 10...\n"
		if [[ -f "$6"_3_mapped.bam ]]
			then echo -e "Output file already existing. Skiping step\n" 
		else
			samtools view -h -F 4 -q 10 -b "$6"_3.bam > "$6"_3_mapped.bam
		fi
		
		echo -e "\nRemoving PCR duplicates...\n"
		if [[ -f "$6"_3_mapped_rmdup.bam ]]
			then echo -e "Output file already existing. Skiping step\n" 
		else
			samtools rmdup -s "$6"_3_mapped.bam "$6"_3_mapped_rmdup.bam
		fi
		
		echo -e "\nStarting convertion in SAM file with samtools..."
		if [[ -f "$6"_3_mapped_rmdup.sam ]]
			then echo -e "\nOutput file already existing. Skiping step\n" 
		else
			samtools view -h "$6"_3_mapped_rmdup.bam > "$6"_3_mapped_rmdup.sam
		fi
		cd ..
		
		mkdir -p step4
		cd step4
		
		echo -e "\nIndexing BAM file..."
		samtools index ../step3/"$6"_3_mapped_rmdup.bam
		echo -e "\nExtractSClip..."
		extractSClip.pl -i ../step3/"$6"_3_mapped_rmdup.bam --ref_genome ../"$3"
		
		cd ..
		
	# Creation d'un bedgraph pour visualisation sur UCSC
	
	#	echo -e "\n\tStarting generating BedGraph..." 
	#	if [[ -f "$6".bedGraph ]]
	#		then echo -e "Output file already existing. Skiping step\n" 
	#	else
	#		genomeCoverageBed -ibam "$6"_3_mapped_rmdup.bam -bg -trackline -trackopts "name="$6" color=250,0,0" > "$6".bedGraph
	#	fi
	#	

		
	# Creation d'un rapport contenant les flagstats
		echo -e "\n\tGenerating a flagstat report"
		
		if [[ -f "$6"_report.txt ]]
			then echo -e "Output file already existing. Skiping step\n" 
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
			
			echo -e "\n\tReference 3 path = "$3"" >> "$6"_report.txt
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
		
			echo -e "\nPrinting the report\n"
			cat "$6"_report.txt
		fi
	
	# rangement dossier
	#	gzip step*/*.fastq
	
		exit 0
	
	else 
		echo -e "\t6 parameters are required \n"
		echo -e "\t1 = Reference 1 Bowtie2 index-file"
		echo -e "\t2 = Reference 2 Bowtie2 index-file"
		echo -e "\t3 = Reference 3 BWA index-file"
		echo -e "\t4 = R1.fastq(.gz)"
		echo -e "\t5 = R2.fastq(.gz)"
		echo -e "\t6 = Output.name\n"
		exit 1
	fi
