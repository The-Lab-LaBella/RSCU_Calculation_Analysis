#!/bin/bash
#SBATCH --time=300:00:00
#SBATCH --partition=Orion
#SBATCH --mem=94GB
#SBATCH --mail-user=bzavalam@uncc.edu
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL

module load hisat2/
module load samtools/
module load anaconda3/2022.10
module load transdecoder/
conda activate metagenomics

#Building an index for reference genomes
# hisat2-build -p 16 hanseniaspora_uvarum.fas hanseniaspora_uvarum_genome

# hisat2-build -p 16 yHMPu5000034959_hanseniaspora_occidentalis_var_occidentalis_170307.fas hanseniaspora_occidentalis_genome


# directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/occid"

# #Loop through the forward reads
# for forward_read in "$directory"/*-R1-trimmed.fastq.gz; do
# 	#Extract the filename without the extension 
# 	filename="${forward_read%-R1-trimmed.fastq.gz}"
# 	#Construct the filename of the corresponding reverse read
# 	reverse_read="${filename}-R2-trimmed.fastq.gz"
	
# 	#Perform operations of the paired-end files
# 	# Perform operations on the paired-end files
#     echo "Forward read: $forward_read"
#     echo "Reverse read: $reverse_read"

#     # Add your desired commands here using "$forward_read" and "$reverse_read"
# 	hisat2 -q -x hanseniaspora_occidentalis_genome -1 $forward_read -2 $reverse_read -S $filename.sam
	
# 	#convert and sort the sam into bam files
# 	samtools view -bS -o $filename.bam $filename.sam
# 	samtools sort -o $filename.bam $filename.bam
	
# done

# uva_directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/uvarum"

# #Loop through the forward reads
# for forward_read in "$uva_directory"/*-R1-trimmed.fastq.gz; do
# 	#Extract the filename without the extension 
# 	filename="${forward_read%-R1-trimmed.fastq.gz}"
# 	#Construct the filename of the corresponding reverse read
# 	reverse_read="${filename}-R2-trimmed.fastq.gz"
	
# 	#Perform operations of the paired-end files
# 	# Perform operations on the paired-end files
#     echo "Forward read: $forward_read"
#     echo "Reverse read: $reverse_read"

#     # Add your desired commands here using "$forward_read" and "$reverse_read"
# 	hisat2 -q -x hanseniaspora_uvarum_genome -1 $forward_read -2 $reverse_read -S $filename.sam
	
# 	#convert and sort the sam into bam files
# 	samtools view -bS -o $filename.bam $filename.sam
# 	samtools sort -o $filename.bam $filename.bam
	
# done


# ##Gene (transcript) Annotation
# directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/occid"

# for alignment_file in "$directory"/*.bam; do

# 	filename="${alignment_file%.bam}"
# 	#annotate transcripts from refrence annotation file
# 	#-C is used for full covergae, end to end by reads
# 	stringtie -o $filename.gtf -G /projects/labella_lab/all_files_y1000_plus/y1000p_gtf_files/yHMPu5000034959_hanseniaspora_occidentalis_var_occidentalis_170307.final.gtf $alignment_file

# done



# uva_directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/uvarum"

# for alignment_file in "$uva_directory"/*.bam; do

# 	filename="${alignment_file%.bam}"
# 	#-C is used for full covergae, end to end by reads
# 	#annotate transcripts from refrence annotation file
# 	stringtie -o $filename.gtf -G /projects/labella_lab/all_files_y1000_plus/y1000p_gtf_files/hanseniaspora_uvarum.final.gtf $alignment_file

# done


##Extract cds from gtf files
#https://github.com/TransDecoder/TransDecoder/wiki
# directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/occid"

# for gtf_file in "$directory"/*.gtf; do
# 	filename="${gtf_file%.gtf}"
# 	gtf_genome_to_cdna_fasta.pl $gtf_file yHMPu5000034959_hanseniaspora_occidentalis_var_occidentalis_170307.fas > $filename.CDS.fasta
# done

# uva_directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/uvarum"
# for gtf_file in "$uva_directory"/*.gtf; do
# 	filename="${gtf_file%.gtf}"
# 	gtf_genome_to_cdna_fasta.pl $gtf_file hanseniaspora_uvarum.fas > $filename.CDS.fasta
# done


##expression count data
# stringtie --merge -p 8 -G /projects/labella_lab/all_files_y1000_plus/y1000p_gtf_files/yHMPu5000034959_hanseniaspora_occidentalis_var_occidentalis_170307.final.gtf -o occid_merged.gtf merged_gtf_occid.txt
# stringtie --merge -p 8 -G /projects/labella_lab/all_files_y1000_plus/y1000p_gtf_files/hanseniaspora_uvarum.final.gtf -o uvarum_merged.gtf merged_gtf_uvarum.txt


##Gene (transcript) Annotation
directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/occid"

for alignment_file in "$directory"/*.bam; do

	filename="${alignment_file%.bam}"
	#annotate transcripts from refrence annotation file
	#-C is used for full covergae, end to end by reads
	stringtie -e -B -G /projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/assembly/occid_merged.gtf $alignment_file -o $filename_abun.gtf

done



uva_directory="/projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/for_alnment/uvarum"

for alignment_file in "$uva_directory"/*.bam; do

	filename="${alignment_file%.bam}"
	#-C is used for full covergae, end to end by reads
	#annotate transcripts from refrence annotation file
	stringtie -e -B -G /projects/labella_lab/all_files_y1000_plus/hanseniaspora_sequencing_data/assembly/uvarum_merged.gtf $alignment_file -o $filename_abun.gtf

done