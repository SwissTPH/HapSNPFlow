######################
# Script for merging all the read processing summaries in order
# to extract summary statistics
#
# 07.07.2023
# monica.golumbeanu@unibas.ch
######################

library(ggplot2)

all_sample_file = "/scicore/home/pothin/golmon00/GitRepos/sara_ampseq/sample_metadata/sample_table.txt"
deplex_marker_summary_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Sara/output_files/SummaryDePlexMarker/"
merge_summary_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Sara/output_files/SummaryMerge/"

all_sample_tab = read.csv(all_sample_file, header = TRUE)

# Read the demultiplexing results
merged_tab_deplex = NULL
for (summary_file_deplex in list.files(deplex_marker_summary_folder, full.names = TRUE)) {
  file_tab = read.table(summary_file_deplex, sep = "\t", header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
  merged_tab_deplex = rbind.data.frame(merged_tab_deplex, file_tab)
}

# Identify which samples were dropped after demultiplexing
for (marker_id in unique(merged_tab_deplex$MarkerID)) {
  dropped_samples = setdiff(unique(all_sample_tab$SampleID), 
                            unique(merged_tab_deplex[which(merged_tab_deplex$MarkerID == marker_id), "SampleID"]))
  print(paste("Samples which were dropped after demultiplexing for ", marker_id,  ":"))
  print(paste(dropped_samples, collapse = ","))
}
   
# Read the merging results
merged_tab_merge = NULL
for (summary_file_merge in list.files(merge_summary_folder, full.names = TRUE)) {
  file_tab = read.table(summary_file_merge, sep = "\t", header = TRUE, stringsAsFactors = FALSE, as.is = TRUE)
  merged_tab_merge = rbind.data.frame(merged_tab_merge, file_tab)
}

# Identify which samples were dropped after merging
all_dropped = NULL
for (marker_id in unique(merged_tab_merge$MarkerID)) {
  dropped_samples = setdiff(unique(all_sample_tab$SampleID), 
                            unique(merged_tab_merge[which(merged_tab_merge$MarkerID == marker_id), "SampleID"]))
  print(paste("Samples which were dropped after merging for ", marker_id,  ":"))
  print(paste(dropped_samples, collapse = " "))
  all_dropped = c(all_dropped, dropped_samples)
}
all_dropped = all_dropped[which(!grepl("NTC", all_dropped))]

# Merge all results for the remaining samples
all_tab = merge(merged_tab_deplex, merged_tab_merge, by = c("SampleID", "BarcodePair", "SampleName", "MarkerID") )

all_tab$proc_deplexed = all_tab$numReadOut/all_tab$numReadIn
all_tab$proc_merged = all_tab$numRead/all_tab$numReadOut

# Print details of missing samples
all_dropped_tab = merged_tab_deplex[which(merged_tab_deplex$SampleID %in% all_dropped),]

# Plot distribution of read counts after demultiplexing
ggplot(all_tab, aes(x = proc_deplexed*100, group = MarkerID)) + 
  geom_histogram(binwidth=1) +
  theme_bw() +
  facet_wrap(~ MarkerID) + xlab ("% reads assigned after demultiplexing") +
  ylab("Number of samples")

# Plot distribution of read counts after merging
ggplot(all_tab, aes(x = proc_merged*100, group = MarkerID)) + 
  geom_histogram(binwidth=1) +
  theme_bw() + 
  facet_wrap(~ MarkerID) + xlab ("% reads assigned after merging") +
  ylab("Number of samples")


