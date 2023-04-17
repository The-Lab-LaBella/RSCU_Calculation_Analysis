#!/bin/bash
#SBATCH --partition=Orion
#SBATCH --time=15:00:00

#This is a simple Bash script just to add commas seperating the columns to convert them into a csv file

for file in *out;
do
	awk '{print $1","$2","$3","$4","$5","$6}' $file >> $file.csv
	sed -i '1i Sequence,Codon,Amino acid,Frequency,Percentage,Codon Count' $file.csv
done
