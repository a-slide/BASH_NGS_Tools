#! /bin/bash
    
    if [[ $# != 3 ]]; then
        prog_name="$(basename "$(test -L "$0" && readlink "$0" || echo "$0")")"
        echo -e "\nUsage: $prog_name <ref.fasta> <seq.bam> <output_prefix>"
        echo -e "\t* fasta needs to be uncompressed"
        echo -e "\t* bam needs to be sorted and indexed\n"
        exit 1
    fi
    
    echo -e "Sort and index files"
    samtools faidx $1
    
    if [[ -f "$3".bcf ]]
    then echo -e "\nBCF file already exists\n" 
    else
        echo -e "\nCompute BCF file\n" 
        samtools mpileup -g -L 1000000 -d 1000000 -o 10 -e 10 -f $1 $2 > "$3".bcf
    fi
    
    echo -e "\nFilter and convert to vcf\n"
    # Remove ref and soft clip lowQual and LowDepth
    bcftools call -c -O v < "$3".bcf |\
    bcftools filter -e "%TYPE='ref'" - |\
    bcftools filter -s "LowQual" -e "%QUAL<30" - |\
    bcftools filter -s "LowDepth" -e "INFO/DP<200" - |\
    bcftools filter -s "NoDepth" -e "INFO/DP<20" - > "$3"_filtered.vcf
    
    #echo -e "\nGenerate stats and graphs\n"
    #bgzip -i -f -c "$3"_filtered.vcf >"$3".vcf.gz
    #tabix -p vcf $3.vcf.gz
    #bcftools stats -f $1 -s - "$3".vcf.gz > "$3"_stat.txt
    #mkdir -p plots/
    #plot-vcfstats -p plots/ "$3"_stat.txt
    
    echo -e "\nCreate consensus fastq sequence\n"
    # Keep all ref but hard clip low depth
    bcftools call -c -O v < "$3".bcf |\
    bcftools filter -e "INFO/DP<20" - |\
    vcfutils.pl vcf2fq > "$3".fastq
    
    exit 0
