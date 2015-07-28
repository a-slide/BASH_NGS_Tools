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

function controlcreation(){
	
	fastq_control_sampler -e 0.005 -n 10000 $(tail -n +3 refseq.txt | awk '{print $2"/"$1".fa 10000"}' | tr "\n" " ")
	echo "Creation of fastq controls (R1 + R2) : Success."
	control=1
}

########################################################################
#                                                                      #
#	MAIN SCRIPT                                                        #
#                                                                      #
########################################################################

# Option parsing

OPTS=$( getopt -o ht:f:r: -l,bowtie:,create_control,preprocessing,quality_control,edit_refseq_list -- "$@" )

eval set -- "$OPTS"

while true; do
	case $1 in
		-h) printhelp
		exit 0;;

		-t) nbt="$2"
		shift 2;;

		--edit_refseq_list) editrefseqlist
		shift;;
		
		--create_control) 
			controlcreation
		shift;;
		
		--preprocessing) preprocessing
		shift;;

		--quality_control) qualitycontrol
		shift;;
		
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

# reading of the options config file : refseq.txt
nbLignes=`wc -l ./refseq.txt | cut -f1 -d" "`
mkdir -p results/
tmpdir1=$(mktemp -d)

# Mapping loop
for i in $(seq 3 $nbLignes); do
	bn=`head -n "$i" ./refseq.txt | tail -n 1 | cut -f1` 			# basename of fasta file
	path=`head -n "$i" ./refseq.txt | tail -n 1 | cut -f2`			# path of fasta file
	indexPath=`head -n "$i" ./refseq.txt | tail -n 1 | cut -f3`		# path of index file
	
	echo -e "\n########################################################"
	echo -e "# Mapping against $bn"
	echo -e "########################################################\n"
	
	# Preparation of the folder
	mkdir -p $bn
	cd $bn
	tmpdir2=$(mktemp -d)
	
	mv ../*R1.fastq ./
	mv ../*R2.fastq ./

	echo -e "### Aligning with Bowtie2 ###\n"      
	bowtie2 --local --"$sens" -X 2000 --reorder -x "$indexPath"/"$bn" -p "$nbt" -1 ./*R1.fastq -2 ./*R2.fastq | samtools view - -hbS | samtools sort - "$tmpdir2"/"$bn"
	
	echo -e "### Processing alignment files ###\n" 
	samtools view -h "$tmpdir2"/"$bn".bam | head -n 1000 | grep '^@' > "$tmpdir2"/header.sam &
	samtools view -h -q 30 "$tmpdir2"/"$bn".bam | cut -f 1 | grep -v '^@' | sort | uniq > "$tmpdir2"/mapped.txt &
	samtools view -h "$tmpdir2"/"$bn".bam | cut -f 1 | grep -v '^@' | sort | uniq > "$tmpdir2"/all.txt &
	samtools view -h -f 64 "$tmpdir2"/"$bn".bam | grep -v '^@'| sort > "$tmpdir2"/R1.sam &
	samtools view -h -f 128 "$tmpdir2"/"$bn".bam | grep -v '^@'| sort > "$tmpdir2"/R2.sam
	wait
	
	comm -13 "$tmpdir2"/mapped.txt "$tmpdir2"/all.txt > "$tmpdir2"/unmapped.txt
	wait
	
	cat "$tmpdir2"/header.sam > $tmpdir1/"$bn".sam
	join -t $'\t' "$tmpdir2"/mapped.txt "$tmpdir2"/R1.sam >> $tmpdir1/"$bn".sam
	join -t $'\t' "$tmpdir2"/mapped.txt "$tmpdir2"/R2.sam >> $tmpdir1/"$bn".sam 
	wait
	
	# if last iteration copy the unmapped sequence to results folder
	if [ "$bn" == $(tail -n 1 ../refseq.txt | cut -f1) ]; then
		echo -e "### Exporting unmapped sequences for blastn ###\n" 
		join -t $'\t' "$tmpdir2"/unmapped.txt "$tmpdir2"/R1.sam | grep -v "^@" >> $tmpdir1/temp.sam
		join -t $'\t' "$tmpdir2"/unmapped.txt "$tmpdir2"/R2.sam | grep -v "^@" >> $tmpdir1/temp.sam
		sort $tmpdir1//temp.sam -o $tmpdir1//temp.sam
		wait
	# else generate fastq for the next iteration
	else 
		echo -e "### Generating fastq for the next iteration ###\n" 
		join -t $'\t' "$tmpdir2"/unmapped.txt "$tmpdir2"/R1.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ../unmapped_"$bn"_R1.fastq &
		join -t $'\t' "$tmpdir2"/unmapped.txt "$tmpdir2"/R2.sam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ../unmapped_"$bn"_R2.fastq
		wait
	fi

	mv "$tmpdir2"/*.txt "$tmpdir2"/"$bn".bam ./
	rm -rf $tmpdir2
	cd ..

done;


echo -e "\n########################################################"
echo -e "# Refinement of mapping with blast"
echo -e "########################################################\n"

cd ./results/

	fa=`head -n 3 ../refseq.txt | tail -n 1 | cut -f1`
	ref=`head -n 3 ../refseq.txt | tail -n 1 |  cut -f2`

	# Blast
	awk '{print ">"$1"\n"$10}' $tmpdir1/temp.sam | blastn -subject "$ref"/"$fa".fa -evalue 10e-20 -dust no -outfmt 6 | awk 'BEGIN{FS="\t"} $7 > 99 {print $1}'| sort | uniq > mapped.txt

	# Remove hits from the last rev and add them to the first ref
	cut -f 1 $tmpdir1/temp.sam | grep -v '^@' | sort | uniq > all.txt
	comm -13 mapped.txt all.txt > unmapped.txt
	
	join -t $'\t' mapped.txt "$tmpdir1"/temp.sam >> "$tmpdir1"/"$fa".sam &
	join -t $'\t' unmapped.txt $tmpdir1/temp.sam > "$tmpdir1"/unmapped.sam
	wait
	echo "`grep '^@' "$tmpdir1"/"$fa".sam; grep -v "^@" "$tmpdir1"/unmapped.sam | sort`" > "$tmpdir1"/unmapped.sam
	rm -f $tmpdir1/temp.sam
	
echo -e "\n########################################################"
echo -e "# Creating output files"
echo -e "########################################################\n"


echo -e "\n## Creating diagramms ##\n\n"
	for i in $tmpdir1/*.sam; do
		echo "$(basename $i .sam);$(sed '/^@.*/d' $i | wc -l)" >> ./freqmapped.csv
	done;
	MAW_graph.R ./freqmapped.csv &> /dev/null

echo -e "\n## Creating bam and bedgraphs for visualisation with igv ##\n\n"
	
for i in $(seq 3 $nbLignes); do
	bn=`head -n "$i" ../refseq.txt | tail -n 1 | cut -f1` # basename of references
	samtools view -hbS "$tmpdir1"/"$bn".sam | samtools sort - "$bn"
	samtools index "$bn".bam &
	genomeCoverageBed -ibam "$bn".bam -bg -trackline -trackopts "name="$bn" color=250,0,0" | awk '$4 > 1'> "$bn".bedGraph
	wait
done;
	
echo -e "\n## Creating a list of overrepresented genomic DNA loci > 5 ##\n\n"
	tail -n +2 "$bn".bedGraph | awk '$4 > 5' | sort -r -V -k 4 > genomic_DNA_loci.txt
	
rm -rf $tmpdir1

echo -e "\n########################################################"
echo -e "# DONE"
echo -e "########################################################\n"








