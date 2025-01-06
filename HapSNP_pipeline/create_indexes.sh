#!/bin/bash
#SBATCH --job-name=create_index
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.o
#SBATCH --time=00:30:00
#SBATCH --qos=30min

########################
# Create indexes for the markers
# 02.05.2024
# monica.golumbeanu@swisstph.ch
#
# sbatch create_indexes.sh 
#######################
module purge
module load bwakit
module load SAMtools

output_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/indexes/"

# Create the output folder
if [ -d "$output_folder" ]; then
    echo "Folder $output_folder already exists."
else
    # Create the folder
    mkdir "$output_folder"
    echo "Folder $output_folder created."
fi

# marker information
marker_info_file="/scicore/home/pothin/golmon00/GitRepos/ampseq_pipeline/analysis_Daniela_Paragon_run1/Marker_information.txt"

# Loop through each line in the file
while IFS= read -r line; do
    # Extract the first, second, and third columns
    column1=$(echo "$line" | awk '{print $1}')
    column2=$(echo "$line" | awk '{print $2}')
    column3=$(echo "$line" | awk '{print $3}')

    # Print the extracted columns (optional)
    # echo "First column: $column1, Second column: $column2, Third column: $column3"
    echo "Indexing $column1"

    index_prefix_file=$output_folder
    reference_genome_file=$output_folder$column1.fasta
    # Perform further processing with the extracted columns 
    # Generate index files for bwa
    # bwa index -p $index_prefix_file"/"$column1 $reference_genome_file 
    # Generate .fai file
    samtools faidx $reference_genome_file
done < "$marker_info_file"
