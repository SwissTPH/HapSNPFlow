#!/bin/bash
#SBATCH --job-name=ampSeq_analysis
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=15G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_callHap_%A.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_callHap_%A.o
#SBATCH --qos=1day

########################
# Submit the SNP and haplotype calling part for all samples together
# 07.07.2023
# monica.golumbeanu@swisstph.ch
#
# sbatch submit_AmpSeqCallHaplotypes.sh /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/ /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/
#######################

module purge
module load R/4.3.2-foss-2023a

input_folder=$1
output_folder=$2

echo "Calling SNPs and haplotypes"
Rscript haplotypR_callSNPs_haplotypes.R $input_folder $output_folder




