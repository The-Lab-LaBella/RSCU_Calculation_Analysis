#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=Orion

for file in *csv;
do
        #Remove first line of the file
        sed -i '1d' $file
        #append the codon and its RSCU average value to a new file
        awk '{print $2","$3}' $file > Edited_$file
        #Include their perspective headers in that file
        sed -i '1i Codon,Mean' Edited_$file
done
