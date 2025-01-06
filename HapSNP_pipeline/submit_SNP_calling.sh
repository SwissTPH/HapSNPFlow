#!/bin/bash
#SBATCH --job-name=SNP_calling
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/SNP_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/SNP_%A_%a.o
#SBATCH --qos=6hours

#############################
# Call SNPs in all samples
#
# monica.golumbeanu@unibas.ch
# 05.10.2024
#
# Usage sbatch --array=1-63 submit_SNP_calling.sh /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/reference/PlasmoDB-62_Pfalciparum3D7_Genome.fasta /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/aligned/ /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/SNPs/
#############################

# Inputs
reference_genome_file=$1
aligned_folder=$2
out_folder=$3

# Load necessary modules
module purge
module load BCFtools

# For testing
# reference_genome_file=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/reference/PlasmoDB-62_Pfalciparum3D7_Genome.fasta
# aligned_folder=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/aligned/
# out_folder=/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/SNPs/

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
echo $selected_aligned_file
# Extract sample name
sample_name=$(basename "${selected_aligned_file%%_sorted.bam}")
echo $sample_name
aligned_sample_file=$selected_aligned_file
sample_vcf_file=$out_folder$sample_name"_SNPs.vcf"
sample_vcf_file_all=$out_folder$sample_name"_all.vcf"
echo $sample_vcf_file

bcftools mpileup -Ou -f $reference_genome_file -a AD,DP -Q 0 -q 0 --min-BQ 0 --max-depth 10000000 $aligned_sample_file | bcftools call -m -A -Ov -o $sample_vcf_file_all

bcftools mpileup -Ou -f $reference_genome_file -a AD,DP --min-BQ 0 -Q 0 -q 0 --max-depth 10000000 $aligned_sample_file | bcftools call -mv -Ov -o $sample_vcf_file
