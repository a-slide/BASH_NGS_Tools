#! /bin/bash

########################################################################
# Paramètres
# 1 = chemin du fichier contenant les numéros PDB des proteines d'interet
# 2 = identifiant du groupe (nom de sortie)
########################################################################

if [[ $# == 2 ]]; then

#Test si paramètre est un fichier ou un dossier
#Si fichier = dl
#Sinon utilisation des pdb dans le dossier.


## TELECHARGEMENT DES FICHIERS PDB #####################################

mkdir -p PBD_"$2"/
cd PBD_"$2"/
echo -e "\nDowloading files from PDB"
for file in `cat ../$1`
do
    echo -e "\t$file"
    wget "http://www.rcsb.org/pdb/files/"$file".pdb"
done

## CREATION DES FICHIERS DSSP ##########################################

mkdir -p ../DSSP_"$2"
echo -e "\nCreating DSSP file from PDB files"
for file in *.pdb
do
echo -e "\t"$file""
	temp=`echo "$file" | sed "s/\(.*\)\.pdb/\1/"` # enregistre le nom du fichier sans .pdb dans une variable temporaire
	dssp "$file" > ../DSSP_"$2"/"$temp".dssp
done

## EXTRACTION DES SEQUENCES ET STRUCTURES ###############################

cd ../DSSP_"$2"/
mkdir -p ../MATRIX_"$2"
rm -f ../MATRIX_"$2"/"$2"_AA_SS.txt

echo -e "\nExtracting sequences and structures from DSSP files"
for file in *.dssp
do
	echo -e "\t$file"
	temp=`echo $file | sed "s/\(.*\)\.dssp/\1/"` # enregistre le nom du fichier sans .dssp dans une variable temporaire
	echo -e "$temp" >> ../MATRIX_"$2"/"$2"_AA_SS.txt
	tail -n +29 "$file" | cut -c 14 | tr "[:lower:]" "[:upper:]" | tr "!" "-" | tr -d "\n" >> ../MATRIX_"$2"/"$2"_AA_SS.txt
	echo -e "" >> ../MATRIX_"$2"/"$2"_AA_SS.txt
	tail -n +29 "$file" | cut -c 17 | tr "[:lower:]" "[:upper:]" | tr "!" "-" | tr -d "\n" | tr "[:space:]" "C" >> ../MATRIX_"$2"/"$2"_AA_SS.txt
	echo -e "" >> ../MATRIX_"$2"/"$2"_AA_SS.txt
done

##  PREPARATION DES MATRICES DE FREQUENCE ##############################

cd ../MATRIX_"$2"/
ls -l

./matrice "$2"_AA_SS.txt
# Creation de 3 matrices de fréquences
# AA_SS.csv		Frequence de chaque SS en fonction de l'AA à la position concernée
# AA+1_SS.csv	Frequence de chaque SS en fonction de l'AA à la position suivante
# AA-1_SS.csv	Frequence de chaque SS en fonction de l'AA à la position precedante


## USAGE ###############################################################

else
	echo -e "Wrong number of arguments."
	echo -e "\tArg1 [Input path of PDB files]"
	echo -e "\tArg2 [Output name]"
	exit 1 
fi
