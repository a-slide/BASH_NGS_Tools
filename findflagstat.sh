#! /bin/zsh

for i in `find . -name "*$1"`
do
		echo "processing file : $i \n"
		samtools flagstat "$i" > "$i".flagstat.txt
done
