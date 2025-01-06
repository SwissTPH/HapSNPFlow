#!/bin/bash
#SBATCH --job-name=demux_aligned
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/demux_aligned_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/demux_aligned_%A_%a.o
#SBATCH --qos=6hours

#############################
# Demultiplex reads based on
# their overlap with regions
#
# monica.golumbeanu@unibas.ch
# 07.09.2024
#
# Usage sbatch --array=1-63 submit_demultiplex_aligned.sh /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/marker_info/ /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/aligned/ /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/demultiplexed_aligned/
#############################

# Inputs
BED_folder=$1
aligned_folder=$2
out_folder=$3

# For testing
BED_folder=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/marker_info/
aligned_folder=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/aligned/
out_folder=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/demultiplexed_aligned/

# Create the output folder if it does not exist
if [ -d "$out_folder" ]; then
    echo "Folder $out_folder already exists."
else
    # Create the folder
    mkdir "$out_folder"
    echo "Folder $out_folder created."
fi

# Retrieve array ID
# IDX=1
IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)

# List all aligned files (bam) and store them in an array
aligned_files=($(find "$aligned_folder" -type f -name '*.sorted.bam' | sort))

# Select the NIDX-th file (IDX is 1-based, so subtract 1 for 0-based index)
selected_aligned_file=${aligned_files[IDX-1]}

# Extract sample name
sample_name="$(basename "${selected_aligned_file%_sorted.bam}")"

# List all refion files (bed) and store them in an array
region_files=($(ls "$BED_folder"/*.bed))

# echo "${region_files[@]}"

# Loop through all the region files
for region_file in "${region_files[@]}"
do
  # Extract marker name
  marker_name="$(basename "${region_file%.*}")"
  echo "Processing $marker_name"

  # Extract the aligned reads that overlap the selected region
  demultiplexed_aligned_name=$out_folder$sample_name"_"$marker_name".bam"
  echo "Intersecting with marker region ..."
  module purge
  module load BEDTools
  bedtools intersect -a $selected_aligned_file -b $region_file -f 0.90 -wa > $demultiplexed_aligned_name

  # Extract the forward and reverse reads corresponding to that region
  demultiplexed_f_reads=$out_folder$sample_name"_"$marker_name"_R1.fastq.gz"
  demultiplexed_r_reads=$out_folder$sample_name"_"$marker_name"_R2.fastq.gz"
  demultiplexed_reads=$out_folder$sample_name"_"$marker_name".fastq.gz"
  echo "Extracting reads ..."
  module purge
  module load SAMtools
  # samtools fastq "$demultiplexed_aligned_name" > "$demultiplexed_reads"

  samtools fastq -1 "$demultiplexed_f_reads" -2 "$demultiplexed_r_reads" -0 /dev/null -s /dev/null -n "$demultiplexed_aligned_name"
  echo "Done demultiplexing"

  # Remove unnecessary files
  # rm $demultiplexed_aligned_name
done
