###########################
# Script using the HaplotypR package to call SNPs and detect 
# haplotypes in AmpSeq data
# 
# created 07.07.2023
# monica.golumbeanu@unibas.ch
##########################

# Load the necessary packages
library("HaplotypR")
library("ShortRead")
library(dplyr)

# Retrieve inputs
args = commandArgs(TRUE)
input_folder = args[1]
output_folder = args[2]
# 
# input_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/"
# output_folder = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/output_files/"

# Load the marker sequences (Forward, Reverse, Reference sequence)
primerFile = file.path(paste0(input_folder, "markerFile.txt"))
markerTab = read.delim(primerFile, stringsAsFactors = F)

# Load the combined file with details on the processed reads
procReads = read.table(paste0(output_folder, "merged_reads_table.txt"), sep="\t", header = TRUE)

# Build the reference sequences
postfix = "_merge"
refSeq = DNAStringSet(markerTab$ReferenceSequence)
names(refSeq) = markerTab$MarkerID
lapply(seq_along(refSeq), function(i){
  writeFasta(refSeq[i], file.path(output_folder, paste(names(refSeq)[i], postfix, ".fasta", sep="")))
})

###############
# Calling SNPs
###############

print("Calling SNPs ...")

# Calculate mismatch rate and call SNPs

# Options
minMMrate = 0.05
minOccGen = 1 # modified because we only have one replicate but it should occur in 2 out of 3!

# For debugging single markers
# markerTab = markerTab %>% filter(MarkerID == "marker_10_pfcarl2" | MarkerID == "marker_1_cpmp1")

# process each marker
snpLst = lapply(markerTab$MarkerID, function(marker){
  # Calculate mismatch rate
  seqErrLst = calculateMismatchFrequencies(as.character(procReads[procReads$MarkerID == marker, "ReadFile"]), 
                                            refSeq[marker], 
                                            method ="pairwiseAlignment", # c("pairwiseAlignment","compareDNAString"), 
                                            minCoverage = 100L) 
  names(seqErrLst) = procReads[procReads$MarkerID == marker, "SampleID"]
  seqErr <- do.call(cbind, lapply(seqErrLst, function(l){
    l[,"MisMatch"]/l[,"Coverage"]
  }))
  write.table(seqErr, file.path(output_folder, sprintf("mismatchRate_rate_%s%s.txt", marker, postfix)), sep="\t", row.names=F)
  
  # Call SNPs
  potSNP = callGenotype(seqErr, minMismatchRate = minMMrate, minReplicate = minOccGen)
  snpRef = unlist(lapply(potSNP, function(snp){
    as.character(subseq(refSeq[marker], start=snp, width=1))
  }))
  
  if(!isEmpty(potSNP)){
    snps = data.frame(Chr = marker, Pos = potSNP, Ref = snpRef, Alt = "N", stringsAsFactors = F)
    write.table(snps, file = file.path(output_folder, sprintf("potentialSNPlist_rate%.0f_occ%i_%s%s.txt",
                                                              minMMrate*100, minOccGen, marker, postfix)),
                row.names = F, col.names = T, sep = "\t", quote = F)
    
    # Plot mismatch rate and SNP calls
    png(file.path(output_folder, sprintf("plotMisMatchRatePerBase_rate%.0f_occ%i_%s%s.png",
                                         minMMrate*100, minOccGen, marker, postfix)),
        width=1500 , height=600, type = "cairo")
    matplot(seqErr, type = "p", pch = 16, cex = 0.4, col = "#00000088", ylim = c(0, 1),
            ylab = "Mismatch Rate", xlab = "Base Position", main = marker, cex.axis = 2, cex.lab = 2)
    abline(v = snps[,"Pos"], lty = 2, col = "grey")
    abline(h = minMMrate, lty = 1, col = "red")
    dev.off()
  } else {
    snps = NULL
  }

  return(snps)
})
names(snpLst) = markerTab$MarkerID

###############
# Calling Haplotypes
###############

print("Calling haplotypes ...")

# call haplotype options
minCov = 3
detectionLimit = 1/100
minOccHap = 2
minCovSample = 25

# Select only markers that had SNPs detected
null_elements = names(snpLst)[sapply(snpLst, is.null)]
markerTab = markerTab %>% filter(!(MarkerID %in% null_elements))

# call final haplotypes
finalTab <- createFinalHaplotypTable(
  outputDir = output_folder, sampleTable = procReads, markerTable = markerTab, referenceSeq = refSeq,
  snpList = snpLst, postfix = postfix, minHaplotypCoverage = minCov, minReplicate = minOccHap, 
  detectability = detectionLimit, minSampleCoverage = minCovSample)

saveRDS(finalTab, file.path(output_folder, "finalHaplotypeTab.RDS"))

