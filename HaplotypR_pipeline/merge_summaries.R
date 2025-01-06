######################
# Script for merging all the read processing summaries; this file 
# will be used for calling SNPs and detecting haplotypes
#
# 26.03.2023
# monica.golumbeanu@unibas.ch
######################

# Example of output folder for debugging:
# out_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/output_files/"

# Retrieve inputs
args = commandArgs(TRUE)
out_folder = args[1]

# Build the paths of the relevant folders
merge_summary_folder = paste0(out_folder, "/SummaryMerge/")
merge_summary_file = paste0(out_folder, "/merged_reads_table.txt")

merged_tab = NULL
for (summary_file in list.files(merge_summary_folder, full.names = TRUE)) {
  file_tab = read.table(summary_file, sep = "\t", header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
  merged_tab = rbind.data.frame(merged_tab, file_tab)
}

write.table(merged_tab, merge_summary_file, sep= "\t", row.names = FALSE, quote = FALSE)
