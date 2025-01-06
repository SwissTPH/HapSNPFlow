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
# For Paragon run 1 there are 63 samples (2 NTCs)
# sbatch --array=1-62 submit_alignment.sh  
#######################
module purge
module load BWA

IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)
# IDX=10

output_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/aligned/"
trimmed_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt/"
index_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/reference/PlasmoDB-62_Pfalciparum3D7_Genome.fasta"
sample_file="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files_old/sample_names.csv"

# Create the output folder
if [ -d "$output_folder" ]; then
    echo "Folder $output_folder already exists."
else
    # Create the folder
    mkdir "$output_folder"
    echo "Folder $output_folder created."
fi

# Select the IDX sample
sample_name=$(sed "${IDX}q;d" $sample_file)
echo $sample_name

# Identify the files
in_file_R1=$trimmed_folder$sample_name"_R1.fastq.gz"
in_file_R2=$trimmed_folder$sample_name"_R2.fastq.gz"
out_file_aligned=$output_folder$sample_name".sam"

echo "$in_file_R1 $in_file_R2 $out_file_aligned"

# Perform alignment
bwa mem $index_folder $in_file_R1 $in_file_R2 > $out_file_aligned

# Convert to a bam file, sort and index
module purge
module load SAMtools
bam_file_aligned=$output_folder$sample_name".bam"
bam_file_sorted=$output_folder$sample_name".sorted.bam"
samtools view -b $out_file_aligned > $bam_file_aligned;
samtools sort $bam_file_aligned > $bam_file_sorted;
samtools index $bam_file_sorted
rm $out_file_aligned

# Alignment stats
raw_read_count=$(zcat $in_file_R2 | wc -l)
raw_read_count=$((raw_read_count / 4))
avg_read_count=$((raw_read_count / 62))
aligned_reads=$(samtools view -c -F 4 $bam_file_aligned)
log_file=$output_folder$sample_name".log"
echo "$sample_name $raw_read_count $avg_read_count $aligned_reads" > $log_file

