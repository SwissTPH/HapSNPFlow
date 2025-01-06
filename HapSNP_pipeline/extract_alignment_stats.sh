#!/bin/bash
#SBATCH --job-name=alignment_stats
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.o
#SBATCH --time=00:30:00
#SBATCH --qos=30min
#
####################
# extract reads mapped and unmapped per sample
#
# created 03.05.2024
# monica.golumbeanu@swisstph.ch
# sbatch --array=1-2060 extract_alignment_stats.sh
####################
module purge
module load bwakit

IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)
# IDX=10

alignment_folder="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_panel_redesign/processed_reads/aligned/"
alignment_files=$(ls $alignment_folder/*.sorted.bam)
output_file="/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_panel_redesign/alignment_stats/file_"$IDX".txt"
# echo $output_file

folders=$(echo "$alignment_files" | sed 's/ /\n/g')
sample_folder=$(echo "$folders" | awk -v N="$IDX" 'NR == N')
# echo $sample_folder

aligned_n=$(samtools stats $sample_folder | grep 'SN' | grep 'reads mapped:' | cut -f3)
unmapped_n=$(samtools stats $sample_folder | grep 'SN' | grep 'reads unmapped:' | cut -f3)
echo $sample_folder" "$aligned_n" "$unmapped_n > $output_file

