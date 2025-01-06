#!/bin/bash
#SBATCH --job-name=align
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/align_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/align_%A_%a.o
#SBATCH --time=00:30:00
#SBATCH --qos=6hours

########################
# Demultiplex reads by marker
# 02.05.2024
# monica.golumbeanu@swisstph.ch
#
# sbatch --array=1-64 submit_alignment.sh 
#######################
module purge
module load bwakit

# IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)
IDX=10
output_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/aligned/"
trimmed_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt/"
index_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/indexes/"
# marker information
marker_info_file="/scicore/home/pothin/golmon00/GitRepos/ampseq_pipeline/analysis_Daniela_Paragon_run1/Marker_information.txt"

# Create the output folder
if [ -d "$output_folder" ]; then
    echo "Folder $output_folder already exists."
else
    # Create the folder
    mkdir "$output_folder"
    echo "Folder $output_folder created."
fi

# identify the samples
dirs=$(find /scicore/projects/openbis/userstore/uni_basel_stph_nsanzabana/ -mindepth 1 -maxdepth 1 -type d -newermt 2022-11-21 ! -newermt 2022-11-22)

# Extract relevant files
folders=$(echo "$dirs" | sed 's/ /\n/g')
sample_folder=$(echo "$folders" | awk -v N="$IDX" 'NR == N')
# echo $folders

in_file_R1=$(ls $sample_folder/*_R1_*.fastq.gz)
in_file_R2=$(ls $sample_folder/*_R2_*.fastq.gz)

# Loop through each line in the file
while IFS= read -r line; do
    # Extract the first, second, and third columns
    column1=$(echo "$line" | awk '{print $1}')
    column2=$(echo "$line" | awk '{print $2}')
    column3=$(echo "$line" | awk '{print $3}')

    # Print the extracted columns (optional)
    echo "First column: $column1, Second column: $column2, Third column: $column3"

    # Perform further processing with the extracted columns
    R1_with_ext=$(basename "$in_file_R1")
    sample_name="${R1_with_ext%.fastq.gz}"
    input_file_R1=$trimmed_folder$sample_name".fastq.gz"
    R2_with_ext=$(basename "$in_file_R2")
    sample_name="${R2_with_ext%.fastq.gz}"
    input_file_R2=$trimmed_folder$sample_name".fastq.gz"
    echo $input_file_R1

    modified_sample_name=$(echo "$sample_name" | sed 's/_R2_//g')
    out_file_aligned=$output_folder$modified_sample_name"_"$column1".sam"
    index_prefix=$index_folder$column1
    reference_file=$index_folder$column1".fasta"
    echo $out_file_aligned
    
    # Perform alignment
    bwa mem $index_prefix $input_file_R1 $input_file_R2 > $out_file_aligned
    
    # Convert to a bam file, sort and index
    bam_file_aligned=$output_folder$modified_sample_name"_"$column1".bam"
    bam_file_sorted=$output_folder$modified_sample_name"_"$column1".sorted.bam"
    samtools view -b $out_file_aligned > $bam_file_aligned;
    samtools sort $bam_file_aligned > $bam_file_sorted;
    samtools index $bam_file_sorted
    rm $out_file_aligned
    
    # Alignment stats
    raw_read_count=$(zcat $input_file_R2 | wc -l)
    raw_read_count=$((raw_read_count / 4))
    avg_read_count=$((raw_read_count / 59))
    aligned_reads=$(samtools view -c -F 4 $bam_file_aligned)
    log_file=$output_folder$modified_sample_name"_"$column1".log"
    echo "$sample_name $raw_read_count $avg_read_count $aligned_reads" > $log_file
done < "$marker_info_file"
