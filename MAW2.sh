#!/bin/bash
export LC_ALL=C

# variables:
nbt=1
sens="sensitive"
control=0

# FUNCTIONS
function printhelp(){
	echo -e "Welcome in the usage of Multiple Alignment Workflow MAW\nCall:\n\t./MAW [OPTIONS] ./R1file.fastq.gz ./R1file.fastq.gz\n\nOPTIONS:\n"
    echo -e "\t-h: provide help, don't need fastq file"
    echo -e "\t--edit_refseq_list: call gedit in order to fill ./refseq.txt"
    echo -e "\t\tFrom the third line and for each line corresponding to one reference genome, the format is:\n\t\t basename of reference genome <tabulation> absolute path of bowtie2 indexes <tabulation> absolute path of fasta file of the reference genome"
    echo -e "\t--preprocessing: call sickle in order to do triming on input data"
    echo -e "\t--quality_control: call FastQC to do a quality control on input data"
    echo -e "\t--nb_thread: to define the number of useful threads for alignements by Bowtie2"
    echo -e "\t--bowtie: to define the sensitivity of alignement by Bowtie2 (accepted values: very-fast, fast, sensitive and very-sensitive)"
    echo -e "\t--create_control: create data set controls with default values"
    echo -e "You must have R"
    echo -e "You must add Bowtie2, Samtools the rscript.R to your PATH"
    echo -e "If required, you might need to add FastQC and sickle to your PATH"
}

function preprocessing(){
	if [ "$(which sickle)" == "" ]; then
		echo "Intall Sickle, or add it to your PATH."
		exit 0
	fi
	
	echo "Pre-processing on-going ..."
	mkdir -p preprocessing
	mv ./R1.fastq ./preprocessing
	mv ./R2.fastq ./preprocessing
	cd ./preprocessing
	
	sickle pe -f R1.fastq -r R2.fastq -t sanger -q 30 -l 100 -o trim_R1.fastq -p trim_R2.fastq -s trim_single.fastq
	
	mv ./trim_R1.fastq ..
	mv ./trim_R2.fastq ..

	echo "Pre-processing done."
}

function qualitycontrol(){
	if [ "which fastqc*" == "" ]; then
		echo "Intall Fastqc, or add it to your PATH."
		exit 0
	fi
	echo "Quality control on-going ..."
	fastqc *.fastq
	echo "Quality control done."
}

function editrefseqlist(){
	if [ ! -f ./refseq.txt ]; then
		touch ./refseq.txt
	fi
	gedit ./refseq.txt
}

#function index(){
	#indexPath=$1
	#filepath=$2
	#ls $filepath
	#if [ "$(ls "$filepath")" != "" ]; then
		#if [ ! -f "$indexPath"*.bt2 ]; then
			#mkdir -p "$indexPath"
			#cd "$indexPath"
			#bowtie2 build "$filepath" -S .
		#fi
	#else
		#echo "Incorrect path for fasta files."
		#exit 0
	#fi
#}

function controlcreation(){
	
	fastq_control_sampler -e 0.001 -n 100000 $(tail -n +3 refseq.txt | awk '{print $2"/"$1".fa"}' | tr "\n" " ")
	echo "Creation of fastq controls (R1 + R2) : Success."
	control=1
}

########################################################################
#                                                                      #
#	MAIN SCRIPT                                                        #
#                                                                      #
########################################################################

# Option parsing

OPTS=$( getopt -o hf:r: -l nb_thread:,bowtie:,create_control,preprocessing,quality_control,edit_refseq_list -- "$@" )

eval set -- "$OPTS"

while true; do
	case $1 in
		-h) printhelp
		exit 0;;

		--edit_refseq_list) editrefseqlist
		shift;;
		
		--create_control) 
			controlcreation
		shift;;
		
		--preprocessing) preprocessing
		shift;;

		--quality_control) qualitycontrol
		shift;;

		--nb_thread) nbt="$2"
		shift 2;;
		
		--bowtie)
			if [ "$2" == 'very-fast' ] || [ "$2" == 'fast' ] || [ "$2" == 'sensitive' ] || [ "$2" == 'very-sensitive' ]; then
				sens="$2"
			else
				echo "bowtie2 wrong option"
				exit 0
			fi
		shift 2;;

		--)shift;break;;
	esac

done
shift $(($OPTIND - 1))

# Verify dependancies
if [ "$(which bowtie2)" == "" ] || [ "$(which samtools)" == "" ] || [ "$(which bedtools)" == "" ]; then
	echo "Intall Bowtie2 or/and Samtools or/and Bedtools. Or add them to your PATH."
	exit 0
fi

# fastq parameters
if [ "$control" != 1 ] ; then
	if [ "$1" == "" ] || [ "$2" == "" ]; then
		echo "You must precise R1.fastq.gz and R2.fastq.gz files."
		exit 1
	else
		echo -e "### Decompressing Fastq files ###\n"
		gunzip -c "$1" > ./R1.fastq
		gunzip -c "$2" > ./R2.fastq
	fi
fi

nbLignes=`wc -l ./refseq.txt | cut -f1 -d" "`
# reading of the options config file : refseq.txt
mkdir -p mapped_reads/

# Mapping loop
for i in $(seq 3 $nbLignes); do
	bn=`head -n "$i" ./refseq.txt | tail -n 1 | cut -f1` 			# basename of fasta file
	path=`head -n "$i" ./refseq.txt | tail -n 1 | cut -f2`			# path of fasta file
	indexPath=`head -n "$i" ./refseq.txt | tail -n 1 | cut -f3`		# path of index file
	
	echo -e "\n########################################################"
	echo -e "\n# Mapping against $bn"
	echo -e "\n########################################################\n\n"
	mkdir -p $bn
	mv ./*R1.fastq $bn
	mv ./*R2.fastq $bn

	cd $bn

	echo -e "### Aligning with Bowtie2 ###\n"      
	bowtie2 --local --"$sens" -X 2000 --reorder -x "$indexPath"/"$bn" -p "$nbt" -1 ./*R1.fastq -2 ./*R2.fastq -S "$bn".sam
	
	samtools view -hS -q 30 "$bn".sam | cut -f 1 | grep -v '^@' | sort | uniq > MR_"$bn".txt
	
	samtools view -hS -f 64 "$bn".sam | sort > R1_"$bn".sam
	samtools view -hS -f 128 "$bn".sam | sort > R2_"$bn".sam
	
	grep '^@' "$bn".sam > ../mapped_reads/MR_"$bn".sam
	join -t $'\t' MR_"$bn".txt R1_"$bn".sam >> ../mapped_reads/MR_"$bn".sam
	join -t $'\t' MR_"$bn".txt R2_"$bn".sam >> ../mapped_reads/MR_"$bn".sam
	
	cut -f 1 "$bn".sam | grep -v '^@' | sort | uniq > ALL_"$bn".txt
	comm -13 MR_"$bn".txt ALL_"$bn".txt > UR_"$bn".txt
	
	join -t $'\t' UR_"$bn".txt R1_"$bn".sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ../UR_"$bn"_R1.fastq
	join -t $'\t' UR_"$bn".txt R2_"$bn".sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ../UR_"$bn"_R2.fastq
	
	cd ..

done;

########################################################################
# Creating output files
########################################################################

echo -e "Output data managing ...\n"

cd ./mapped_reads/

echo -e "\n########################################################"
echo -e "\n# Creating a list of overrepresented genomic DNA loci"
echo -e "\n########################################################\n\n"

cut -f 3,4 MR_"$bn".sam | sort | rev | sed 's/^...\(.*\)/\1/g' | rev | uniq -f 1 -c > ./genomic_DNA_loci.txt

join -t $'\t' ../"$bn"/UR_"$bn".txt ../"$bn"/R1_"$bn".sam >> temp.sam
join -t $'\t' ../"$bn"/UR_"$bn".txt ../"$bn"/R2_"$bn".sam >> temp.sam

sort temp.sam -o temp.sam

echo -e "\n########################################################"
echo -e "\n# BLAST : Unmapped"
echo -e "\n########################################################\n\n"

fa=`head -n 3 ../refseq.txt | tail -n 1 | cut -f1`
ref=`head -n 3 ../refseq.txt | tail -n 1 |  cut -f2`

tail -n +4 temp.sam | cut -f1,10 | sed 's/^/>/g' | sed 's/\t/\n/g' | blastn -subject "$ref"/"$fa".fa -evalue 10e-20 -dust no -outfmt 6 | awk 'BEGIN{FS="\t"} $7 > 99 {print $1}'| sort | uniq > mapped.txt
cut -f 1 temp.sam | grep -v '^@' | sort | uniq > all.txt

comm -13 mapped.txt all.txt > unmapped.txt

join -t $'\t' mapped.txt temp.sam >> MR_$fa.sam
join -t $'\t' unmapped.txt temp.sam > unmapped.sam

echo "`grep '^@' MR_"$bn".sam; grep -v "^@" unmapped.sam | sort`" > temp.sam && mv temp.sam unmapped.sam

echo -e "\n########################################################"
echo -e "\n# Creating diagramms"
echo -e "\n########################################################\n\n"

for i in *.sam; do
	echo "$(basename $i .sam | sed 's/MR_//g');$(sed '/^@.*/d' $i | wc -l)" >> ./freqmapped.csv
done;

MAW_graph.R ./freqmapped.csv
echo -e "done.\n"

for i in *.sam; do
	samtools view -hbS "$i" | samtools sort - "$(basename $i .sam)"
	samtools index "$(basename $i .sam).bam"
done;



#if [ "$control" == 1 ] ; then
	#for i in $(seq 3 $nbLignes); do
		#bn=`head -n "$i" ./refseq.txt | tail -n 1 | cut -f1`
		#grep -c "$bn" 
#fi

##script variant calling
#java -jar GenomeAnalysisTK.jar -T UnifiedGenotyper -L chr$i -I input_file.bam -glm BOTH -R "genome de ref" -o output_file.vcf





