##########################
# Look at alignment logs and print coverage
#
# monica.golumbeanu@swisstph.ch
# 20.08.2024
##########################

library("Rsamtools")
library(ShortRead)
library(dplyr)

align_dir = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/aligned/"
demux_dir = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/demultiplexed_aligned/"

# Sample information
sample_info = read.csv("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/sample_table.txt")
# Marker information
marker_info = read.csv("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/markerFile.txt", sep = "\t")

# Loop through all the log files
result_table = NULL
for (sample_name in unique(sample_info$SampleName)) {
  aligned_bam_file = list.files(align_dir, pattern = paste0(sample_name, "\\.bam$"), full.names = TRUE)

  count_data = tryCatch({
      countBam(aligned_bam_file)$records
  },
  error = function(e) {
      return(0)  # Return NULL or any other appropriate action in case of an error
  },
  warning = function(w) {
      message("Warning encountered: ", w)
      # You can handle warnings specifically if needed
      return(countBam(aligned_bam_file)$records)  # Continue with the operation or handle it
  }
  )

  # Select the demultiplexed marker files
  if(count_data > 0) {
      demux_list = list.files(path = demux_dir, pattern = paste0(sample_name, "\\.sorted\\.bam_", ".*\\.bam$"), full.names = TRUE)
      print(sample_name)
      if (length(demux_list) > 0) {
          for (marker_file in demux_list) {
              marker_name = strsplit(marker_file, ".bam")[[1]][2]
              count_data_marker = countBam(marker_file)
              row_stats = c(sample_name, marker_name, floor(count_data/2), floor(count_data_marker$records/2))
              result_table = rbind.data.frame(result_table, row_stats)

          }
      }
  }

}

# Change column names of the final table
colnames(result_table) = c("sample_name", "marker", "total_aligned_reads", "demultiplexed_reads")
result_table$marker = paste0("marker", result_table$marker)
result_table$marker = gsub("\\.", "_", result_table$marker)

# Load results from HaplotypR
hap_results = read.csv("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/plots/plots_read_counts/Read_count_summaries.csv")
hap_results = hap_results %>% select(SampleID, Total, MarkerID, Merged) %>% rename(sample_name = SampleID,
                                                                                               marker = MarkerID)
# Merge all results
all_summaries = merge(result_table, hap_results, by = c("sample_name", "marker"), all = TRUE)
all_summaries = all_summaries %>% rename(total_demultiplexed_HaplotypR = Total,
                         total_merged_HaplotypR = Merged,
                         total_aligned_demultiplexed = demultiplexed_reads)
write.csv(all_summaries, "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/plots/plots_read_counts/summaries_comparing_approaches.csv")



