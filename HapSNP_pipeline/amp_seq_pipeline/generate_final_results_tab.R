##################################
# Script for generating the final results 
#
# created 10.07.2023
# monica.golumbeanu@unibas.ch
##################################

final_hap_tab = readRDS("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Sara/output_files2/finalHaplotypeTab.RDS")

output_folder = "/scicore/home/pothin/golmon00/genotyping/sara_ampseq/"

write.table(final_hap_tab$csp, file.path(output_folder, "raw_results_TabHaplotype_csp.txt"), sep="\t", row.names=F, quote=F)
write.table(final_hap_tab$cpmp, file.path(output_folder, "raw_results_TabHaplotype_cpmp.txt"), sep="\t", row.names=F, quote=F)
write.table(final_hap_tab$`ama1-D3`, file.path(output_folder, "raw_results_TabHaplotype_ama1-D3.txt"), sep="\t", row.names=F, quote=F)
write.table(final_hap_tab$`ama1-D2`, file.path(output_folder, "raw_results_TabHaplotype_ama1-D2.txt"), sep="\t", row.names=F, quote=F)
write.table(final_hap_tab$cpp, file.path(output_folder, "raw_results_TabHaplotype_cpp.txt"), sep="\t", row.names=F, quote=F)
