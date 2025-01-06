#!/bin/bash
#SBATCH --job-name=ampSeq_analysis
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.o
#SBATCH --time=00:30:00
#SBATCH --qos=30min

########################
# Submit ampSeq analysis for all samples
# 26.03.2024
# monica.golumbeanu@swisstph.ch
#
# sbatch --array=1-131 submit_AmpSeqPreprocess.sh /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/ /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/output_files/
#######################

module load R/4.1.0-foss-2018b

IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)
input_folder=$1
output_folder=$2

# Create the folders for different steps of the pipeline
# Folder with results for demultiplexing by marker
demultiplex_marker_folder=$output_folder"/dePlexMarker/"
if [ -d "$demultiplex_marker_folder" ]; then
    echo "Folder $demultiplex_marker_folder already exists."
else
    # Create the folder
    mkdir "$demultiplex_marker_folder"
    echo "Folder $demultiplex_marker_folder created."
fi

# Folder with results for merging reads
merge_folder=$output_folder"/processedReads/"
if [ -d "$merge_folder" ]; then
    echo "Folder $merge_folder already exists."
else
    # Create the folder
    mkdir "$merge_folder"
    echo "Folder $merge_folder created."
fi

echo "Demultiplexing reads by marker ..."
Rscript haplotypR_demultiplex_marker_sample.R $input_folder $output_folder $IDX

echo "Merge forward and reverse reads"
Rscript haplotypR_merge_reads_sample.R $input_folder $output_folder $IDX




