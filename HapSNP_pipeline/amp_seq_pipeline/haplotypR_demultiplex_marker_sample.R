###########################
# Script using the HaplotypR package to demultiplex and merge
# paired end AmpSeq reads
# 
# created 07.07.2023
# monica.golumbeanu@unibas.ch
##########################

# Load the necessary packages
library("HaplotypR")
library("ShortRead")

print(Sys.time())
start_time = Sys.time()

# Retrieve inputs
args = commandArgs(TRUE)
input_folder = args[1]
output_folder = args[2]
idx = as.integer(args[3])

# Input and output folders:
# input_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Sara/input_files/"
# output_folder = "/scicore/home/pothin/golmon00/genotyping/sara_ampseq/"
# idx = 10
# output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Sara/output_files/"

##################
# LOADING SAMPLE
##################
## Construct the file paths of the input files necessary for the workflow and 
# check if the files exist

print("Locating input files...")

# File with the marker sequences (Forward, Reverse, Reference sequence)
primerFile = file.path(paste0(input_folder, "markerFile.txt"))
if(!file.exists(primerFile)) {
  stop("Primer file not found. Make sure that a file called markerFile.txt exists
       in the input folder")
}

# Sample file containing the details of the demultiplexed samples
dePlexSample = read.csv(file.path(paste0(input_folder, "sample_table.txt")))
# Select only one sample for testing
dePlexSample = dePlexSample[idx, ]

################
# Demultiplexing
################

print(paste("Demultiplexing", dePlexSample$SampleID, " by marker ..."))

# Demultiplexing by marker 
# OUTPUT: 
# dePlexMarker/ folder contains for each sample and marker the corresponding reads (.fastq.gz)
# summary_dePlexMarker/SummaryDemultiplexMarker_idx.txt contains for each sample, marker and barcode pair the number of reads and the paths to the files containing the reads

# Output subdirectory to store the reads after demultiplexing by marker
outDeplexMarkerDir = file.path(output_folder, "dePlexMarker", idx)
if (dir.exists(outDeplexMarkerDir)) {
  print("Removing already existing folder for marker demultiplexing.")
  unlink(outDeplexMarkerDir, recursive = TRUE)
}
dir.create(outDeplexMarkerDir)

# Output subdirectory to store the summaries after demultiplexing by marker
outDeplexMarkerSummaryDir = file.path(output_folder, "SummaryDePlexMarker")
dir.create(outDeplexMarkerSummaryDir, showWarnings = FALSE)

# Process each marker
markerTab = read.delim(primerFile, stringsAsFactors = F)
dePlexMarker = demultiplexByMarker(dePlexSample, markerTab,
                                   outDeplexMarkerDir,
                                   trimFilenameExt = "R1_001_MM_1\\.fastq.gz$",
                                   max.mismatch = 2)
dePlexMarker = dePlexMarker[!is.na(dePlexMarker$FileR1),]

# Save summary table
write.table(dePlexMarker, file.path(outDeplexMarkerSummaryDir, paste0("SummaryDemultiplexMarker_", idx, ".txt")),
            sep="\t", row.names = F, quote = FALSE)

