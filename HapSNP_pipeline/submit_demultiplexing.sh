#!/bin/bash
#SBATCH --job-name=demultiplex
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.o
#SBATCH --qos=6hours

########################
# Demultiplex reads by marker
# 02.05.2024
# monica.golumbeanu@swisstph.ch
#
# sbatch --array=1-63 submit_demultiplexing.sh /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt/ /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/barcodes_F.fasta /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/barcodes_R.fasta /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/demultiplexedMarker/
# sbatch --array=1-131 submit_demultiplexing.sh /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/cutadapt/ /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/barcodes_F.fasta /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/barcodes_R.fasta /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/demultiplexedMarker/
#######################
module purge
module load cutadapt

IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)

# folder cutadapt/ with the reads where the illumina adapter was removed
input_folder=$1
# file with forward primer information
primer_file_F=$2
# file with reverse primer information
primer_file_R=$3
# folder where the demultiplexed reads will be saved
output_folder=$4

# Create the output folder
if [ -d "$output_folder" ]; then
    echo "Folder $output_folder already exists."
else
    # Create the folder
    mkdir "$output_folder"
    echo "Folder $output_folder created."
fi

# Identify the sample
# List all *R1*.fastq.gz files and store them in an array
R1_files=($(find "$input_folder" -type f -name '*L001_R1_001*.fastq.gz' | sort))

# Select the NIDX-th file (IDX is 1-based, so subtract 1 for 0-based index)
selected_R1_file=${R1_files[IDX-1]}
filename_R1=$(basename -- "$selected_R1_file")
out_name_R1="${filename_R1%.*.*}"
# echo $out_name_R1

# Construct the pair
selected_R2_file="${selected_R1_file/L001_R1_001/L001_R2_001}"
filename_R2=$(basename -- "$selected_R2_file")
out_name_R2="${filename_R2%.*.*}"
# echo $out_name_R2

# Saving logs
log_file=$output_folder$out_name_R1".log"

# Demultiplex
cutadapt -e 1 -g ^file:$primer_file_F -G ^file:$primer_file_R -o $output_folder$out_name_R1-{name}.fastq.gz -p $output_folder$out_name_R2-{name}.fastq.gz $selected_R1_file $selected_R2_file

