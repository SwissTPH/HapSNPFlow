############################ 
# Test DeepSNV
# 
# 25.08.2024
# monica.golumbeanu@unibas.ch
############################

library(deepSNV)
library(Rsamtools)

regions = data.frame(chr = "Pf3D7_11_v3", start = 1294624, stop = 1294798)
Malmix = deepSNV(test = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/aligned/1M_1_sorted.bam",
                  control = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/aligned/PosC10_sorted.bam",
                  regions = regions, q = 10)
                  
# Extract chromosome names to make sure we use the right nomenclature
# # Specify the path to your BAM file
# bam_file <- "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/analysis_Paragon_run1/processed/aligned/1M_1_sorted.bam"
# # Load BAM file header
# bam_header <- scanBamHeader(bam_file)
# # Extract chromosome names
# chromosome_names <- names(bam_header[[1]]$targets)
# # Print chromosome names
# print(chromosome_names)

show(Malmix)
control(Malmix)[100:110,]
test(Malmix)[100:110,]
plot(Malmix)
SNVs <- summary(Malmix, sig.level=0.05, adjust.method="BH")
head(SNVs)
nrow(SNVs)



