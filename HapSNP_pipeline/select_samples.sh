#!/bin/bash

# this script screens the openbis folder and extracts the sample names
# 18.05.2024
# bash select_samples.sh /scicore/projects/openbis/userstore/uni_basel_stph_nsanzabana/ 2022-11-21 /scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_panel_redesign/HaplotypeR/input_files/

openbis_folder=$1
date_created=$2
output_folder=$3

output_file1=$output_folder"/samples_openbis.txt"
output_file2=$output_folder"/sample_table_raw.txt"

echo "SampleName,ForwardReadsFile,ReverseReadsFile,NbForwardReads,NbReverseReads"> $output_file1
echo "SampleID,SampleName,BarcodePair,ReadNumbers,FileR1,FileR2"> $output_file2

echo "Selecting samples created on $date_created"

# Find the directories where your FASTQ files are stored
directories=$(find $openbis_folder -mindepth 1 -maxdepth 1 -type d -newermt "$date_created" ! -newermt "$date_created + 1 day")

# Loop through each directory
for dir in $directories
do
    echo "Processing directory: $dir"
    
    metadata_file="$dir/*metadata.tsv"
    # echo $metadata_file
    # Extract the sample name
    
    sample_name=$(grep '^EXTERNAL_SAMPLE_NAME' $metadata_file | awk '{print $2}')
    echo $sample_name
      
    # Find the forward and the reverse reads
    
    # fastq_R="$dir/*_R2_001*.fastq.gz"
    # fastq_F="$dir/*_R1_001*.fastq.gz"
    fastq_F=$(find "$dir" -type f -name '*_R1_001*.fastq.gz' | head -n 1)
    fastq_R=$(find "$dir" -type f -name '*_R2_001*.fastq.gz' | head -n 1)
    # echo $fastq_R
    
    # Extract the forward and reverse reads
    nb_reads_F=$(zcat $fastq_F | echo $((`wc -l`/4)))
    nb_reads_R=$(zcat $fastq_R | echo $((`wc -l`/4)))
    
    echo $file,$date_created,$sample_name,$fastq_F,$fastq_R,$nb_reads_F,$nb_reads_R
    echo $file,$date_created,$sample_name,$fastq_F,$fastq_R,$nb_reads_F,$nb_reads_R>>$output_file1
    
    echo $sample_name,$sample_name,,$nb_reads_F,$fastq_F,$fastq_R>>$output_file2
done
