########################
# Compute the forward and reverse barcodes
#
# monica.golumbeanu@swisstph.ch
########################

marker_info = read.table("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/markerFile.txt",
                         header = TRUE)
# marker_info = read.table("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/markerFile.txt", 
#                          header = TRUE)
# Write to file
barcode_F_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/barcodes_F.fasta"
barcode_R_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/HaplotypR/input_files/barcodes_R.fasta"
# barcode_F_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/barcodes_F.fasta"
# barcode_R_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_classical_AmpSeq/input_files/barcodes_R.fasta"


if(file.exists(barcode_F_file))
  file.remove(barcode_F_file)
if (file.exists(barcode_R_file))
  file.remove(barcode_R_file)

for (i in 1:nrow(marker_info)){
  write.table(paste0(">", marker_info$MarkerID[i]), barcode_F_file,
              col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  write.table(marker_info$Forward[i], barcode_F_file,
              col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  write.table(paste0(">", marker_info$MarkerID[i]), barcode_R_file,
              col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
  write.table(marker_info$Reverse[i], barcode_R_file,
              col.names = FALSE, row.names = FALSE, quote = FALSE, append = TRUE)
}