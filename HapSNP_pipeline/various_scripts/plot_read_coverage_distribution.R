########################
# Compute histogram of reads per haplotype
# 
# 27.08.2024
# monica.golumbeanu@unibas.ch
########################
library(ggplot2)

results_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Annina/results/finalTabHaplotype_ama1-D3.txt"
haplotype_results = read.table(results_file, sep = "\t", header = TRUE)
haplotype_results = haplotype_results %>% filter(!Haplotype %in% c("Noise", "Singelton", "Chimera", "Indels"))

ggplot(haplotype_results, aes(x = Reads)) +
  geom_histogram(color = "black", fill = "grey", binwidth = 200) +
  labs(title = "", x = "Number of reads", y = "Frequency") +
  theme_bw() 
