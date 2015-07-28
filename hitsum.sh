#! /bin/zsh

for i in `find . -name "*$1"`
	do
		echo "processing file : $i \n"
		echo "Summary of hits \n" > "$i"_hit.txt
		
	for j in {1..22} X Y M
		do
			echo "chromosome $j: " >> "$i"_hit.txt
			cat "$i" | egrep -c "chr"$j"\s" >> "$i"_hit.txt
		done
				
	done
