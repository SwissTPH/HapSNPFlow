#########################################
# create the sample file for the pipeline
#
# monica.golumbeanu@swisstph.ch
# 26.03.2024
##########################################

# Load the details of the samples (folders and the read counts and IDs)
sample_details = read.csv('/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/samples_openbis.txt', 
                          header = FALSE, sep = " ")
# Initialize output table
sample_tab = NULL

# Loop through all the samples and extract the information needed
for (i in 1:nrow(sample_details)) {
  sample_folder = sample_details$V1[i]
  folder_content = list.files(sample_folder, pattern="fastq.gz", full.names = TRUE)
  SampleID = sample_details$V2[i]
  SampleName = substr(SampleID, 1, str_locate(SampleID, "_") - 1)
  BarcodePair = ""
  ReadNumbers = sample_details$V4[i]
  FileR1 = folder_content[1]
  FileR2 = folder_content[2]
  sample_row = cbind(SampleID, SampleName, BarcodePair, ReadNumbers, FileR1, FileR2)
  sample_tab = rbind.data.frame(sample_tab, sample_row)
}

# Write the table to a file
write.csv(sample_tab, "~/GitRepos/sara_ampseq/sample_metadata/sample_table.txt", quote = FALSE, row.names = FALSE)
