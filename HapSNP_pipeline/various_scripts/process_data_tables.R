###############################
# Extract parts of tables
# monica.golumbeanu@swisstph.ch
# 09.07.2023
###############################

tab_markers = read.csv("/scicore/home/mynaco16/GROUP/Teaching/ParagonAmpseq/analysis_Paragon_run3/input/AmpliconPositions/amplicon_positions.csv", sep="\t")
markers_file = "/scicore/home/mynaco16/GROUP/Teaching/ParagonAmpseq/analysis_Paragon_run3/markers"
write_csv(as.data.frame(tab_markers$Amplicon), markers_file, col_names = FALSE)

samples_names = read.csv(samples_file, header = FALSE)
samples_file = "/scicore/home/mynaco16/GROUP/Teaching/ParagonAmpseq/analysis_Paragon_run3/samples"

samples_names_new = gsub("_R[12]", "", samples_names$V1)
write_csv(as.data.frame(samples_names_new), samples_file, col_names = FALSE)

marker_info = read.csv("/scicore/home/mynaco16/GROUP/Teaching/ParagonAmpseq/analysis_Paragon_run3/input/AmpliconPositions/amplicon_positions_old.csv", sep="\t")
amplicon_file = "/scicore/home/mynaco16/GROUP/Teaching/ParagonAmpseq/analysis_Paragon_run3/input/AmpliconPositions/amplicon_positions.csv"
write.table(marker_info, amplicon_file, sep="\t", quote = FALSE, row.names = FALSE)

