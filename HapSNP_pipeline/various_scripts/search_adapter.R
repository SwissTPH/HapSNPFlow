#############################
# Manual search of adapters in a read
#
# monica.golumbeanu@swisstph.ch
# 23.08.2024
#############################

library(dplyr)
library(stringr)
library(ShortRead)

search_adapter = function(adapter_df, read_seq, forward_strand) {
  matching_rows = NULL
  if (forward_strand == "F") {
    matching_rows = adapter_df %>% filter(str_detect(read_seq, Forward))
  } 
  
  if (forward_strand == "R")  {
    matching_rows = adapter_df %>% filter(str_detect(read_seq, Reverse))
  }
  
  if (forward_strand == "Ref")  {
    matching_rows = adapter_df %>% filter(str_detect(read_seq, ReferenceSequence))
  }
  
  return(matching_rows)
}

marker_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/markerFile.txt"
marker_info = read.csv(marker_file, sep= "\t")

# Specify the path to your .fastq.gz file
fastq_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt_demux/BSSE_QGF_223555_000000000_KMM6N_1_CAll2_ACAGGCGC_GTGCGATA_S12_L001_R2_001_MM_1-unknown.fastq.gz"
# Open the FASTQ file for reading
fq = readFastq(fastq_file)

# Loop through each read
all_detected = NULL
for (i in seq_along(fq)) {
  # Extract each read
  read = fq[i]
  
  detected_adapters = search_adapter(marker_info, as.character(sread(read)), "Ref") 
  all_detected = rbind(all_detected, detected_adapters)
}

