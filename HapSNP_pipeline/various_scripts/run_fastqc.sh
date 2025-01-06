#!/bin/bash
#SBATCH --job-name=fastqc
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=10G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_fastqc.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_fastqc.o
#SBATCH --qos=6hours

#########################
# script which runs fastqc on all the .fastq.gz files from a specified list 
# of directories based on date
# sbatch run_fastqc.sh 2022-11-21 /scicore/home/pothin/golmon00/genotyping/Daniela_analyses/Daniela_Paragon_run1/fastqc/
########################

module purge
module load FastQC

# Check if the dateand output folder argument are provided
if [ $# -ne 2 ]; then
    echo "Usage: $0 <date> <output folder>"
    exit 1
fi

# Date to use when selecting folders with results from the sequencing facility
search_date="$1"
# Folder where to save all results
output_folder=$2

# Define the directories where your FASTQ files are stored
# Replace "/path/to/your/directory" with the actual path to your directories
# Store the directories created on the specified date in a variable
directories=$(find /scicore/projects/openbis/userstore/uni_basel_stph_nsanzabana/ -mindepth 1 -maxdepth 1 -type d -newermt "$search_date" ! -newermt "$search_date + 1 day")

# Create the output directory if it doesn't exist
mkdir -p "$output_folder"

# Loop through each directory
for dir in $directories
do
    echo "Processing files in directory: $dir"
    
    # Run FastQC on all FASTQ files in the directory and save results to the output directory
    fastqc "$dir"/*.fastq.gz -o "$output_folder"
done

