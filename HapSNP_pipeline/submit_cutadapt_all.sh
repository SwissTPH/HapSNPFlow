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

# IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)
IDX=10
output_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt_no_dimers/"
output_dimer_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt_dimers/"

# Create the output folder
if [ -d "$output_folder" ]; then
    echo "Folder $output_folder already exists."
else
    # Create the folder
    mkdir "$output_folder"
    echo "Folder $output_folder created."
fi

# Create the output dimer folder
if [ -d "$output_dimer_folder" ]; then
    echo "Folder $output_dimer_folder already exists."
else
    # Create the folder
    mkdir "$output_dimer_folder"
    echo "Folder $output_dimer_folder created."
fi

# identify the samples
dirs=$(find /scicore/projects/openbis/userstore/uni_basel_stph_nsanzabana/ -mindepth 1 -maxdepth 1 -type d -newermt "$date_files" ! -newermt "$date_files + 1 day")

# Extract relevant files
folders=$(echo "$dirs" | sed 's/ /\n/g')
sample_folder=$(echo "$folders" | awk -v N="$IDX" 'NR == N')

in_file_R1=$(ls $sample_folder/*_R1_*.fastq.gz)
in_file_R2=$(ls $sample_folder/*_R2_*.fastq.gz)

# Print the total number of reads
read_count=$(zcat "$in_file_R1" | grep -c "^@")
echo "Total $read_count raw reads."

out_file_R1=$output_folder/$(basename "$in_file_R1")
out_file_R2=$output_folder/$(basename "$in_file_R2")

out_file_dimer_R1=$output_dimer_folder/$(basename "$in_file_R1")
out_file_dimer_R2=$output_dimer_folder/$(basename "$in_file_R2")

log_file=$out_file_R1".log"
cutadapt_json=$out_file_R1".json"

# Remove all adapter dimers
# cutadapt \
#     --action=trim \
#     -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
#     -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
#     -e 0 \
#     --no-indels \
#     --minimum-length 100 \
#     -o $out_file_dimer_R1 \
#     -p $out_file_dimer_R2 \
#     --cores 1 \
#     --json=$cutadapt_json \
#     --compression-level=1 \
#     --quiet \
#     $in_file_R1 \
#     $in_file_R2 > /dev/null
    
# file with forward primer information
primer_file_F="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/barcodes_F.fasta"
# file with reverse primer information
primer_file_R="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/barcodes_R.fasta"

output_folder_demux="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt_demux/"

filename_R1=$(basename -- "$out_file_dimer_R1")
out_name_R1="${filename_R1%.*.*}"
filename_R2=$(basename -- "$out_file_dimer_R2")
out_name_R2="${filename_R2%.*.*}"
cutadapt_json2=$output_folder_demux$out_name_R1".json"

echo $out_name_R1
echo $out_file_dimer_R1

# Print the total number of reads
read_count=$(zcat "$out_file_dimer_R1" | grep -c "^@")
echo "Total $read_count reads without Illumina adapters"

# Demultiplex
cutadapt \
    --action=trim \
    -g ^file:${primer_file_F} \
    -G file:${primer_file_R} \
    --no-indels \
    --trim-n -q 10 \
    --minimum-length 100 \
    -e 0 \
    -o $output_folder_demux$out_name_R1-{name}.fastq.gz \
    -p $output_folder_demux$out_name_R2-{name}.fastq.gz \
    --untrimmed-output $output_folder_demux$out_name_R1-unknown_R1.fastq.gz \
    --untrimmed-paired-output $output_folder_demux$out_name_R2-unknown_R2.fastq.gz \
    --compression-level=1 \
    --quiet \
    $out_file_dimer_R1 \
    $out_file_dimer_R2 > /dev/null

    # --pair-adapters \
    # --discard-untrimmed \
    # -e 0 \
    # --json=${cutadapt_json2} \
    # --untrimmed-output ${trimmed_demuxed_unknown_fastqs}/${sample_id}_unknown_R1.fastq.gz \
    # --untrimmed-paired-output ${trimmed_demuxed_unknown_fastqs}/${sample_id}_unknown_R2.fastq.gz \
total_count=0
while read -r amplicon_name; do
    fastq_file="$output_folder_demux$out_name_R1-${amplicon_name}.fastq.gz"
    if [[ -f $fastq_file ]]; then
        read_count=$(zgrep -c ^@ $fastq_file) || zgrep_exit_status=$?
        zgrep_exit_status=${zgrep_exit_status:-0}
        if [[ $zgrep_exit_status -eq 1 ]]; then
            read_count=0
        elif [[ $zgrep_exit_status -ne 0 ]]; then
            echo "Error: could not calculate amplicons for $fastq_file"
            exit $zgrep_exit_status
        fi
    else
        read_count=0
    fi

    total_count=$((total_count + read_count))
    echo -e "${amplicon_name}\t${read_count}"
    
    # unset the status variable
    unset zgrep_exit_status

done < <(grep '^>' $primer_file_F | sed 's/^>//')
echo $total_count


