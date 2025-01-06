#############################
# Test with DADA2
#
#############################

library(dada2)

path = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/demultiplexed_aligned/" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
#sample.names <- gsub("GFB-8581_000000000-L3VNF_1_0170D212_", "", basename(fnFs))
#sample.names <- gsub("L001_", "", sample.names)
sample.names <- gsub("_R1.fastq.gz", "", basename(fnFs))
sample.names <- gsub(".sorted.bam", "", sample.names)

plotQualityProfile(fnRs[980], aggregate = TRUE)

sample.names.selected = sample.names[1:100]
selected_fnFs = fnFs[1:100]
selected_fnRs = fnRs[1:100]

filtFs <- file.path(path, "filtered", paste0(sample.names.selected, "_filt_R1.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names.selected, "_filt_R2.fastq.gz"))
names(filtFs) <- sample.names.selected
names(filtRs) <- sample.names.selected

out <- filterAndTrim(selected_fnFs, filtFs, selected_fnRs, filtRs, truncLen = c(240, 200),
                     maxN = 0, maxEE = c(5, 5), truncQ = 2, rm.phix=TRUE, #maxEE=c(2, 2)
                     compress = TRUE, multithread = TRUE)
head(out)

filtFs_2 = list.files(file.path(path, "filtered"), pattern = "_filt_R1.fastq.gz")
filtRs_2 = list.files(file.path(path, "filtered"), pattern = "_filt_R2.fastq.gz")

errF <- learnErrors(filtFs_2, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[4]])

seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

