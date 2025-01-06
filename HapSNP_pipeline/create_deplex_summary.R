########################
#
# create summary of demultiplexing by marker
#
# 15.06.2024
# monica.golumbeanu@swisstph.ch
##########################

# Retrieve inputs
args = commandArgs(TRUE)
# sample_file = args[1]
# marker_file = args[2]
# output_folder = args[3]
idx = as.integer(args[1])

# sample_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/sample_table.txt"
# marker_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/markerFile.txt"
# demultiplex_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/demultiplexedMarker/"
# cutadapt_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/cutadapt/"
# output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/"

sample_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/sample_table.txt"
marker_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/markerFile.txt"
demultiplex_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/demultiplexedMarker/"
cutadapt_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/cutadapt/"
output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/processed_reads/"


output_summaries = paste0(output_folder, "SummaryDePlexMarker/")
# idx = 30

# Create output folder
dir.create(output_summaries, showWarnings = FALSE)

# Read in sample and marker information
sample_tab = read.csv(sample_file)
marker_tab = read.csv(marker_file, sep = "\t")

SampleID = sample_tab[idx, "SampleID"]
# SampleID = stringr::str_replace_all(SampleID, "\\.", "_")
SampleID = stringr::str_replace_all(SampleID, "\\.", "")
SampleID = stringr::str_replace_all(SampleID, "_", "")
SampleName = sample_tab[idx, "SampleName"]
# SampleName = stringr::str_replace_all(SampleName, "\\.", "_")
SampleName = stringr::str_replace_all(SampleName, "\\.", "")
SampleName = stringr::str_replace_all(SampleName, "_", "")
SampleR1 = basename(sample_tab[idx, "FileR1"])
SampleR2 = sample_tab[idx, "FileR2"]
BarcodePair = ""

summary_tab = NULL
for (i in 1:nrow(marker_tab)) {
  MarkerID = marker_tab[i, "MarkerID"]
  numReadIn = numReadOut = 0
  
  demultiplexed_files = list.files(demultiplex_folder, full.names = TRUE)
  cutadapt_files = list.files(cutadapt_folder, full.names = TRUE)
  
  # Create the regular expression pattern and search for the files
  pattern_R1 = paste0(SampleID, "_.*", "_R1_", ".*", MarkerID)
  pattern_R2 = paste0(SampleID, "_.*", "_R2_", ".*", MarkerID)
  FileR1 = demultiplexed_files[grep(pattern_R1, demultiplexed_files)]
  FileR2 = demultiplexed_files[grep(pattern_R2, demultiplexed_files)]
  
  # Files before demultiplexing
  pattern1 = paste0(SampleID, "_.*", "_R1_", ".*", ".fastq.gz$" )
  FileR1_AllMarkers = cutadapt_files[grep(pattern1, cutadapt_files)]
  
  # Counting reads
  # Construct the shell commands
  command1 = paste("echo $(zcat", FileR1_AllMarkers, "| wc -l) / 4 | bc")
  command2 = paste("echo $(zcat", FileR1, "| wc -l) / 4 | bc")
  
  # Run the shell command and capture the output
  numReadIn = system(command1, intern = TRUE)
  numReadOut = system(command2, intern = TRUE)
  
  # Concatenate results
  summary_tab = rbind(summary_tab, cbind(SampleID, BarcodePair, SampleName, 
                                         MarkerID, numReadIn, numReadOut,  
                                         FileR1, FileR2))
}
write.table(summary_tab,
            paste0(output_summaries, "SummaryDemultiplexMarker_", idx, ".txt"),
            quote = FALSE, sep="\t", row.names = FALSE)
