########################
# extract read counts for the different steps of the ParagonAmpSeq pipeline
#
# monica.golumbeanu@swisstph.ch
# 17.07.2024
########################
library(stringr)

extract_read_nb = function(input_folder, file_pattern, tab_col_names) {
  # Extract first the raw reads numbers
  file_list = list.files(input_folder, full.names = TRUE, pattern = file_pattern)
  
  read_df = NULL
  for (read_file in file_list) {
    if (grepl(".fastq.gz", file_pattern))
      nb_reads = as.numeric(system(paste("zcat", read_file, "| wc -l | awk '{print $1 / 4}'"), intern = TRUE))
    else
      nb_reads = as.numeric(system(paste("cat", read_file, "| wc -l | awk '{print $1 / 4}'"), intern = TRUE))
    
    sample_name = gsub(file_pattern, "", basename(read_file))
    read_df = rbind.data.frame(read_df, cbind.data.frame(sample_name, nb_reads))
  }
  
  colnames(read_df) = tab_col_names
  return(read_df)
}

extract_read_aligned_nb = function(input_folder, tab_col_names) {
  # Extract first the raw reads numbers
  file_list = list.files(input_folder, full.names = TRUE, pattern = ".sam")
  
  read_df = NULL
  for (read_file in file_list) {
    nb_reads = as.numeric(system(paste("awk '$1 !~ /^@/ && $5 > 0' ", read_file, " | wc -l"), intern = TRUE))
    sample_name = gsub(".sam", "", basename(read_file))
    read_df = rbind.data.frame(read_df, cbind.data.frame(sample_name, nb_reads))
  }
  
  colnames(read_df) = tab_col_names
  return(read_df)
}

extract_read_marker_nb = function(input_folder, file_pattern, tab_col_names) {
  # Extract first the raw reads numbers
  file_list = list.files(input_folder, full.names = TRUE, pattern = file_pattern)
  
  read_df = NULL
  for (read_file in file_list) {
    if (grepl(".fastq.gz", file_pattern))
      nb_reads = as.numeric(system(paste("zcat", read_file, "| wc -l | awk '{print $1 / 4}'"), intern = TRUE))
    else
      nb_reads = as.numeric(system(paste("cat", read_file, "| wc -l | awk '{print $1 / 4}'"), intern = TRUE))
    sample_name = str_sub(basename(read_file), 1, str_locate(basename(read_file), "_marker_")[1] - 1)
    marker_name = str_sub(basename(read_file), str_locate(basename(read_file), "_marker_")[1] + 1, str_locate(basename(read_file), file_pattern)[1] - 1)
    read_df = rbind.data.frame(read_df, cbind.data.frame(sample_name, marker_name, nb_reads))
  }
  
  colnames(read_df) = tab_col_names
  return(read_df)
}

  print("Extracting raw read numbers")
  raw_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/input/raw/"
  # raw_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run3/input/raw/"
  raw_df = extract_read_nb(raw_path, "_R1.fastq.gz", c("sample_name", "total"))

  print("Extracting the number of reads after adapter removal with cutadapt")
  # cutadapt_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run3/processed/cutadapt/"
  cutadapt_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/cutadapt/"
  cutadapt_df = extract_read_nb(cutadapt_path, "_R1.fastq.gz", c("sample_name", "cutadapt"))
  final_df = merge(raw_df, cutadapt_df, by = "sample_name")

  print("Extract the number of reads after merging reads")
  # merged_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run3/processed/fusedReads/"
  merged_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/fusedReads/"
  merged_df = extract_read_nb(merged_path, "_fuseReads.fastq", c("sample_name", "merged"))
  final_df = merge(final_df, merged_df, by = "sample_name")

  print("Extracting the number of reads after aligning the reads")
  # aligned_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run3/processed/aligned/"
  aligned_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/aligned/"
  aligned_df = extract_read_aligned_nb(aligned_path, c("sample_name", "aligned"))
  final_df = merge(final_df, aligned_df, by = "sample_name")
  
  print("Extracting the number of reads after splitting by marker")
  # split_marker_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run3/processed/splitByMarker/"
  split_marker_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/splitByMarker/"
  split_marker_df = extract_read_marker_nb(split_marker_path, ".fastq", c("sample_name", "marker_name", "demultiplexed"))
  final_df = merge(final_df, split_marker_df, by = "sample_name")
  
  print("Extracting the number of reads after removing primers and introns")
  # amplicon_marker_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run3/processed/extractedAmplicons/"
  amplicon_marker_path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/extractedAmplicons/"
  amplicon_df = extract_read_marker_nb(amplicon_marker_path, ".fastq.gz", c("sample_name", "marker_name", "amplicons"))
  final_df = merge(final_df, amplicon_df, by = c("sample_name", "marker_name"))
  
  # write.csv(final_df, "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run3/read_stats_pipeline.csv", 
  #           quote = FALSE, row.names = FALSE)
  
  write.csv(final_df, "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/read_stats_pipeline.csv", 
            quote = FALSE, row.names = FALSE)
  