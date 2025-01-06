#########################
# remove the Pool and one NTC file from sample table
# moncia.golumbeanu@unibas.ch
# 08.09.2024
#########################

sample_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files_old/samples_openbis.txt"
sample_tab = read.table(sample_file, skip = 1, sep = ",")

sample_tab = sample_tab %>% filter(!grepl("Pool", V3))
sample_tab = sample_tab[-13, ]

# Write the sample_names to a one column file
write.table(sample_tab$V3, "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files_old/sample_names.csv",
          quote = FALSE, row.names = FALSE, col.names = FALSE)