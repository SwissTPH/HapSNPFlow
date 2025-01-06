#!/bin/bash
#SBATCH --job-name=cut_adapters
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.o
#SBATCH --time=00:30:00
#SBATCH --qos=30min

########################
# Cut illumina adapters for all samples
# 01.05.2024
# monica.golumbeanu@swisstph.ch
#
# sbatch --array=1-131 submit_cutadapt.sh 2023-05-08 /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/cutadapt
#######################
module purge
module load cutadapt

# Retrieve the date of files and output folder
date_files=$1
output_folder=$2

IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)
# output_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt/"

# Create the output folder
if [ -d "$output_folder" ]; then
    echo "Folder $output_folder already exists."
else
    # Create the folder
    mkdir "$output_folder"
    echo "Folder $output_folder created."
fi

# identify the samples
dirs=$(find /scicore/projects/openbis/userstore/uni_basel_stph_nsanzabana/ -mindepth 1 -maxdepth 1 -type d -newermt "$date_files" ! -newermt "$date_files + 1 day")

# Extract relevant files
folders=$(echo "$dirs" | sed 's/ /\n/g')
sample_folder=$(echo "$folders" | awk -v N="$IDX" 'NR == N')

in_file_R1=$(ls $sample_folder/*_R1_*.fastq.gz)
in_file_R2=$(ls $sample_folder/*_R2_*.fastq.gz)

out_file_R1=$output_folder/$(basename "$in_file_R1")
out_file_R2=$output_folder/$(basename "$in_file_R2")

log_file_R1=$out_file_R1".log"
log_file_R2=$out_file_R2".log"

# run cutadapt
cutadapt -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -o $out_file_R1 $in_file_R1 > $log_file_R1
cutadapt -a AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT -o $out_file_R2 $in_file_R2 > $log_file_R2

