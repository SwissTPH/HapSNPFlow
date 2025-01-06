######################
# Script for merging all the read processing summaries; this file 
# will be used for calling SNPs and detecting haplotypes
#
# 26.03.2023
# monica.golumbeanu@unibas.ch
######################
library(dplyr)

merge_summary_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/SummaryMerge/"
merge_summary_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/merged_reads_table.txt"

# merge_summary_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/SummaryMerge/"
# merge_summary_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/merged_reads_table.txt"

merged_tab = NULL
for (summary_file in list.files(merge_summary_folder, full.names = TRUE)) {
  file_tab = read.table(summary_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
  merged_tab = rbind.data.frame(merged_tab, file_tab)
}

merged_tab = merged_tab %>% filter(numRead > 200)

write.table(merged_tab, merge_summary_file, sep= "\t", row.names = FALSE, quote = FALSE)