#!/bin/bash
#SBATCH --time=72:00:00 
#SBATCH --partition=Orion

#this extracts each sequences from the fasta file and then uses the emboss program to create the codon usage table $1 the reference $2 is the fasta to be parsed

#MUST HAVE CAI conda loaded - must have fasta file all on one line

module load emboss

#ls *.gz | xargs -n1 tar -xzf

#gunzip *gz

for file in *cds;
do
	grep ">" $file | awk '{print $1}' | sed 's/>//' > $file.seq_name 
	
	cat $file.seq_name | while read line
	do
		grep -A 1 "$line" $file > $line.seq.txt

		#do the cusp analysis here
		cusp $line.seq.txt $line.codonTable.out
		
		#you can remove the seq file
		rm $line.seq.txt

		#create one big file for the codon Table data
		#we would need to edit the codon table to remove the lines starting with #
		sed -i "/^#/d" $line.codonTable.out

		#add the sequence name to the first column
		sed -i "s/^/$line /g" $line.codonTable.out

		#add this to a file of ALL the codon tables
		cat $line.codonTable.out >> $file.all_codonTable.out

		rm $line.codonTable.out		
	done
	
done

#cp *all_codonTable.out cds_codonTables
#This is a simple Bash script just to add commas seperating the columns to convert them into a csv file

for file in *out;
do
	awk '{print $1","$2","$3","$4","$5","$6}' $file >> $file.csv
	sed -i '1i Sequence,Codon,Amino acid,Frequency,Percentage,Codon Count' $file.csv
done
