#!/bin/bash
#SBATCH --job-name=summary_demultiplex
#SBATCH --account=mynaco16
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=5G
#SBATCH --output=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.LOG
#SBATCH --error=/scicore/home/pothin/golmon00/genotyping/JOB_OUT/ampSeq_%A_%a.o
#SBATCH --qos=30min

########################
# Demultiplex reads by marker
# 02.05.2024
# monica.golumbeanu@swisstph.ch
#
# sbatch --array=1-63 submit_summary_demultiplexing.sh 

#######################
module purge
module load R

IDX=$(expr ${SLURM_ARRAY_TASK_ID} - 0)

Rscript create_deplex_summary.R $IDX
