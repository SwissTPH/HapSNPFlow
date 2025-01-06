##############################
# Script to create summary of read counts per sample
# 
# monica.golumbeanu@unibas.ch
##############################
library(dplyr)
library(tidyr)
library(ggplot2)

# deplex_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_panel_redesign/HaplotypeR/output_files/SummaryDePlexMarker/"
# merge_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_panel_redesign/HaplotypeR/output_files/merged_reads_table.txt"
# output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_panel_redesign/plots/plots_read_counts/"

# deplex_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/SummaryDePlexMarker/"
# merge_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/merged_reads_table.txt"
# output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/plots/plots_read_counts/"

deplex_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/SummaryDePlexMarker/"
merge_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/merged_reads_table.txt"
output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/plots/plots_read_counts/"


# Create folder where all plots will be stored
dir.create(output_folder, showWarnings = FALSE)

# Create table with demultiplexing summaries
deplex_tab = NULL
for (file_deplex in list.files(deplex_folder, full.names = TRUE)) {
  deplex_sample = read.table(file_deplex, sep = "\t", header = TRUE)
  if(length(colnames(deplex_sample)) == 8) {
    deplex_tab = rbind.data.frame(deplex_tab, deplex_sample)
  }
}
deplex_tab$FileR1 = NULL
deplex_tab$FileR2 = NULL
deplex_tab$BarcodePair = NULL

# Read in the merging results
merge_tab = read.table(merge_file, header = TRUE, sep = "\t")
merge_tab$ReadFile = NULL
merge_tab$BarcodePair = NULL

# Build the final table
final_counts_tab = full_join(deplex_tab, merge_tab, by = c("SampleID", "SampleName", "MarkerID"))
final_counts_tab = final_counts_tab %>% rename(Total = numReadIn,
                                               Demultiplexed = numReadOut,
                                               Merged = numRead)
# Transform into percentages
# final_counts_tab$Demultiplexed = final_counts_tab$Demultiplexed/final_counts_tab$Total
# final_counts_tab$Merged = final_counts_tab$Merged/final_counts_tab$Total
# final_counts_tab$Total = final_counts_tab$Total/final_counts_tab$Total
unique_markers = length(unique(final_counts_tab$MarkerID))

# compute total demultiplexed and total merged
a = as.data.frame(final_counts_tab) %>% 
                      group_by(SampleID, SampleName, Total) %>% 
                      summarise(Total_demultiplexed = sum(Demultiplexed, na.rm = TRUE),
                                Total_merged = sum(Merged, na.rm = TRUE)) %>% ungroup()
a_long = a %>% pivot_longer(-c(SampleID, SampleName), 
                                           names_to = "processing_step", 
                                           values_to = "read_count")
# Build plot with stats per step
# Calculate the percentages
a_long = a_long %>% 
          group_by(SampleID, SampleName) %>% 
          mutate(Percentage = (read_count / read_count[1]) * 100)

# Create table with read numbers for each marker
final_counts_tab_long = final_counts_tab %>% pivot_longer(-c(SampleID, SampleName, MarkerID), 
                                names_to = "processing_step", 
                                values_to = "read_count")
final_counts_tab_long$processing_step = factor(final_counts_tab_long$processing_step, 
                                          levels = c("Total", "Demultiplexed", "Merged"))

# Merge and save counts to file
total_counts_tab = merge(final_counts_tab, a, by = c("SampleID", "SampleName", "Total"))
tab_file = paste0(output_folder, "Read_count_summaries.csv")
write.csv(total_counts_tab, tab_file)
  
# Plot the marker overview for each sample
final_counts_tab_plot = final_counts_tab_long %>% filter(processing_step == "Merged")
p = ggplot(final_counts_tab_plot, aes(x = MarkerID, y = read_count)) +
    geom_boxplot() +
    theme_minimal() +
    theme(legend.position = "none") +
    theme(strip.text.x = element_text(size = 11), 
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  xlab("Marker") + ylab("Read count")
plot_file = paste0(output_folder, "markers_oveview.pdf")
ggsave(plot_file, plot = p,  width = 12, height = 6)

p = ggplot(final_counts_tab_plot, aes(x = MarkerID, y = read_count)) +
  geom_boxplot() +
  theme_minimal() +
  theme(legend.position = "none") +
  theme(strip.text.x = element_text(size = 11), 
        axis.title=element_text(size=10),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
  ylim(c(0,25000))+
  xlab("Marker") + ylab("Read count")
plot_file = paste0(output_folder, "markers_oveview_zoom.pdf")
ggsave(plot_file, plot = p,  width = 12, height = 6)
  
# Plot the reads overview for each sample
for (sample in unique(final_counts_tab_long$SampleName)){
  sample_tab = final_counts_tab_long %>% filter(SampleName == sample)
  sample_overview = a_long %>% filter(SampleName == sample)
  
  sample_overview = sample_overview %>% 
    group_by(SampleID, SampleName) %>% 
    mutate(Average_cov = read_count[1]/unique_markers)
  
  # Create the barplot with the overview for the total reads per sample
  p1 = ggplot(sample_overview, aes(x = processing_step, y = read_count, fill = processing_step)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = paste0(round(Percentage, 2), "%")), 
              vjust = 2, 
              color = "black", 
              size = 3.5) +
    facet_wrap(~SampleID) +
    scale_fill_manual(values=c('#999999','#E69F00', "#2ca25f")) +
    labs(title = "",
         x = "Processing step",
         y = "Count") +
    theme_minimal() +
    theme(legend.position = "none")
  
  # Create the barplot with the overview for the reads per marker
  sample_tab = sample_tab %>% filter(processing_step != "Total")
  sample_tab = merge(sample_tab, unique(sample_overview %>% select(SampleID, SampleName, Average_cov)), 
                     by = c("SampleID", "SampleName"))
  sample_tab$MarkerID = gsub("marker_", "", sample_tab$MarkerID)
  p2 = ggplot(sample_tab, aes(x = MarkerID, y = read_count, fill = as.factor(processing_step))) +
    geom_bar(stat="identity", position = "dodge") +
    geom_hline(aes(yintercept = Average_cov), lty = 2) +
    scale_fill_manual(values=c('#E69F00', "#2ca25f")) +
    facet_wrap(~SampleID) + labs(fill = "Read counts") +
    theme_minimal(base_size = 12) + 
    theme(strip.text.x = element_text(size = 11), 
          axis.title=element_text(size=10),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    xlab("Marker") + ylab("Read count") + theme(legend.position="bottom") 
  
  # Save plot
  plot_file_1 = paste0(output_folder, sample, "_oveview.pdf")
  ggsave(plot_file_1, plot = p1,  width = 10, height = 6)
  plot_file_2 = paste0(output_folder, sample, "_markers.pdf")
  ggsave(plot_file_2, plot = p2,  width = 20, height = 12)
}




