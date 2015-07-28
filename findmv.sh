#! /bin/zsh

for file in `find . -name "*$1"`
do
		cp -rf "$file" "$2" 
done
