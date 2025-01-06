##########################
# Process AmpSeq results - calcuate MOI and haplotype frequencies across replicates
#
# created 12.07.2023
# monica.golumbeanu@unibas.ch
##########################

library(dplyr)
library(stringr)

process_raw_data = function(hap_tab) {
  
  # Ensure the table is a data frame and remove the empty column
  hap_tab = as.data.frame(hap_tab)
  hap_tab$NA. = NULL
  
  # Remove the rows that do not correpsond to a haplotype
  hap_tab = hap_tab %>% dplyr::filter(!(Haplotype %in% c("Noise", "Singelton", "Chimera", "Indels")))
  # Order the rows by sample name
  hap_tab = hap_tab %>% dplyr::arrange(SampleName)
  # Identify the replicate number 
  hap_tab$Replicate = str_sub(hap_tab$SampleID, -1)
  # Calculate MOI for each sample
  hap_tab = hap_tab %>% 
    dplyr::group_by(SampleID) %>% 
    dplyr::mutate(#match(SampleID, unique(SampleID))
                  MOI = length(unique(Haplotype)))
  # Calculate the "frequency" of each haplotype for each replicate
  # To do so, for each replicate, we divide the number of reads of each haplotype 
  # to the total number of reads of the replicate
  hap_tab = hap_tab %>%
    dplyr::group_by(SampleID) %>%
    dplyr::mutate(HapFreqRep = Reads/sum(Reads))
  # We calculate the average frequency for each haplotype across the replicates
  hap_tab = hap_tab %>%
    dplyr::group_by(SampleName, Haplotype) %>%
    dplyr::mutate(HapFreq = mean(HapFreqRep), HapOcc = n())
  
  return(hap_tab)
}

process_marker = function(marker_name, in_folder_name, out_folder_name) {
  
  marker_hap_tab = read.table(paste0(in_folder_name, "raw_results_TabHaplotype_",
                                     marker_name, ".txt"), 
                              stringsAsFactors = FALSE, header = TRUE)
  
  csp_processed = process_raw_data(marker_hap_tab)
  
  write.table(csp_processed, paste0(out_folder_name, "final_results_TabHaplotype_", 
                                    marker_name, ".txt"),
              col.names = TRUE, row.names = FALSE, quote = FALSE)
}

# Process data for all markers
haplotypR_results = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/"
final_results = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/final_results/"
process_marker("csp", haplotypR_results, final_results)
process_marker("cpp",  haplotypR_results, final_results)
process_marker("cpmp",  haplotypR_results, final_results)
process_marker("ama1-D3",  haplotypR_results, final_results)
process_marker("ama1-D2",  haplotypR_results, final_results)
