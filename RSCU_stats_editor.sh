#!/bin/bash
#SBATCH --time=01:00:00
#SBATCH --partition=Orion

for file in *csv;
do
        sed -i '1d' $file
        awk '{print $2","$3}' $file > Edited_$file
        sed -i '1i Codon,Mean' Edited_$file
done
