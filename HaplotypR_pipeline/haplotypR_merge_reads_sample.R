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

# input_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Sara/input_files/"
# output_folder = "/scicore/home/pothin/golmon00/genotyping/sara_ampseq/"
# idx = 10

##################
# LOADING SAMPLE
##################
## Construct the file paths of the input files necessary for the workflow and 
# check if the files exist

print("Locating input files...")

# Load the marker sequences (Forward, Reverse, Reference sequence)
primerFile = file.path(paste0(input_folder, "markerFile.txt"))
markerTab = read.delim(primerFile, stringsAsFactors = F)

# Load the details of the demultiplexed markers
dePlexMarker_file = file.path(paste0(output_folder, "SummaryDePlexMarker/SummaryDemultiplexMarker_", idx, ".txt"))
dePlexMarker = read.table(dePlexMarker_file, sep="\t", header=TRUE)

################
# Merging reads
################
print("Merging reads ...")

## Create output subdirectory for storing the merged reads
outProcFilesDir = file.path(output_folder, "processedReads/")

# Output subdirectory to store the summaries after merging the reads
outMergeSummaryDir = file.path(output_folder, "SummaryMerge")
dir.create(outMergeSummaryDir, showWarnings = FALSE)

postfix = "_merge"
refSeq = DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) = markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(output_folder, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

procReads = mergeAmpliconReads(as.character(dePlexMarker$FileR1), as.character(dePlexMarker$FileR2), 
                               outProcFilesDir, method = "vsearch")

procReads = cbind(dePlexMarker[,c("SampleID", "SampleName","BarcodePair", "MarkerID")], procReads)
# procReads = procReads[-which(procReads$SampleName == "NTC"),]  
procReads = procReads[procReads$numRead >= 190,]
write.table(procReads, file.path(outMergeSummaryDir, paste("processedReadSummary_merge_", idx, ".txt")), sep="\t", row.names=F, quote=F)


