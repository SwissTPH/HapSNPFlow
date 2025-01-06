############################
# Rename the files with adapter cut reads from Paragon run 1 
# using their sample names
#
# monica.golumbeanu@unibas.ch
# 08.09.2024
############################

# Load smaple information
sample_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files_old/samples_openbis.txt"
sample_tab = read.table(sample_file, skip = 1, sep = ",")

# Rename files
cutadapt_f = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt/"
cut_reads_list = list.files(cutadapt_f, pattern = "\\.fastq\\.gz$", full.names = TRUE)

# Remove the pool
cut_reads_list = cut_reads_list[which(!grepl("Undetermined", cut_reads_list))]

# Rename the files
for (cut_f in cut_reads_list) {
  file_basename = basename(cut_f)
  
  # Forward or reverse
  if(grepl("L001_R2_001", cut_f)) {
    strand = "R2"
  } else {
    strand = "R1"
  }
  
  cut_sample = sample_tab %>%
                filter(grepl(file_basename, V4) | grepl(file_basename, V5)) %>%
                select(V3)
  if (nrow(cut_sample) > 1) {
    stop(paste("Same sample name for more than one .fastq.gz pair", cut_sample))
  } else {
    new_file_name = paste0(dirname(cut_f), "/", cut_sample, "_", strand, ".fastq.gz")
    file.rename(from = cut_f, to = new_file_name)
  }
}
