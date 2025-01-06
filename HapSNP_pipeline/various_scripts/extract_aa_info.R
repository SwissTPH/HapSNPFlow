library(VariantAnnotation)
library(Biostrings)
library(rtracklayer)

# Function to transform negative strand exons to positive strand
transform_to_positive = function(exons, genome_seq) {
    # Process each exon
    new_exons = lapply(1:length(exons), function(i) {
        exon = exons[i]
        seq = getSeq(genome_seq, exon)

        if (as.character(strand(exon)) == "-") {
            # Reverse-complement the sequence
            print("negative exon")
            seq <- reverseComplement(seq)

            # Update strand information
            strand(exon) <- "+"
        }

        # Return updated exon (sequence is optional)
        return(exon)
    })

    # Combine updated exons
    return(do.call(c, new_exons))
}

# Get all vcf files
all_vcfs = list.files("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/SNPs/", pattern = ".sorted.bam_SNPs.vcf", full.names = TRUE)
# all_vcfs = list.files("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/SNPs/", pattern = ".sorted.bam_all.vcf", full.names = TRUE)
print("vcf read")
# Read reference genome sequence
fasta_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/reference/PlasmoDB-62_Pfalciparum3D7_Genome.fasta"
reference = readDNAStringSet(fasta_file)
# Extracts the chromosome names so that they are the same as in the vcf
names(reference) = sub("^(\\S+).*", "\\1", names(reference))

# Read annotation file
annotation_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/ParagonAmpseq/reference/PlasmoDB-62_Pfalciparum3D7.gtf"
# Import GTF file
print("import annotation")
gtf = import(annotation_file)
print("annotation imported")

# Filter exons only and transform all to be on the positive strand
exons = gtf[gtf$type == "CDS"]

# transformed_exons = transform_to_positive(exons, reference)
# save(transformed_exons, file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/transformed_exons.Rdata")
load("/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/transformed_exons.Rdata")
# Initialize final data frame with all the SNP and aminoacid information
aa_changes = NULL

print("starting loop")
# Loop through all VCF files
for (vcf_file in all_vcfs) {
    print(vcf_file)
    ### Read VCF file
    # vcf_file = "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/processed_reads/SNPs/1M.1.sorted.bam.vcf"
    vcf = readVcf(vcf_file, "Pf3D7")  # Assuming 'Pf3D7' is your genome assembly
    sample_name = stringr::str_remove_all(basename(vcf_file), ".sorted.bam_all.vcf")
    sample_name = stringr::str_remove_all(basename(vcf_file), ".sorted.bam_SNPs.vcf")

    chromosome_names = seqnames(rowRanges(vcf))

    ### Extract SNP information
    snp_positions = start(rowRanges(vcf))  # Extract positions
    ref_bases = ref(vcf)  # Extract reference bases from VCF object
    alt_bases = alt(vcf)  # Extract alternate bases from VCF object

    # Extract frequencies and quality scores
    qualities <- qual(vcf)  # Quality scores
    depths <- info(vcf)$DP
    allele_depth <- geno(vcf)$AD

    # Step 4: For each variant, extract codons and calculate amino acid change
    for (i in seq_along(snp_positions)) {
        # Get the SNP position, reference, and alternate alleles
        snp_pos = snp_positions[i]
        ref_base = as.character(ref_bases[i])
        if(nchar(as.character(ref_bases[i])) == 1) {
            for (j in 1:(length(alt_bases[[i]]))) {
                alt_base = as.character(alt_bases[[i]][j])  # Take the first alt allele (if multiple)
                if(nchar(alt_base) == 1) {
                    # Determine the chromosome from the VCF
                    chrom = as.character(chromosome_names[i])

                    # Check if the position is within the bounds of the reference sequence
                    if (snp_pos > length(reference[[chrom]])) {
                        warning(paste("SNP position", snp_pos, "is out of bounds for chromosome", chrom))
                        next
                    }

                    # Get the quality score and coverages
                    quality = ifelse(length(qualities) > 0, qualities[i], NA)  # Quality score
                    total_coverage = depths[i]# sum(allele_depth[[i]])  # Total depth
                    ref_coverage = allele_depth[[i]][1]
                    alt_coverage = allele_depth[[i]][j+1]

                    # Find the exon containing the SNP
                    snp_gr = GRanges(seqnames = chrom, ranges = IRanges(start = snp_pos, end = snp_pos))
                    overlapping_exon = transformed_exons[subjectHits(findOverlaps(snp_gr, transformed_exons))]

                    if (length(overlapping_exon) > 1) {
                        print(stop)
                    }

                    # If the SNP is in an exon
                    if (length(overlapping_exon) != 0) {
                        # Rarely, there is one exon on the forward strand and one exon on the reverse strand overlapping
                        for (k in 1:length(overlapping_exon)) {
                            # Extract sequence from the overlapping exon
                            seq = getSeq(reference, overlapping_exon)[k]
                            exon_seq = as.character(seq)

                            # Determine position within codon
                            exon_start = start(overlapping_exon[k])
                            exon_end = end(overlapping_exon[k])
                            relative_pos = snp_pos - exon_start + 1
                            codon_pos = (relative_pos - 1) %% 3 + 1
                            codon_start = relative_pos - (codon_pos - 1)

                            # Extract the codon and apply the SNP
                            codon = substr(exon_seq, codon_start, codon_start + 2)
                            mutated_codon = codon
                            substr(mutated_codon, codon_pos, codon_pos) = alt_base

                            # Translate codons to amino acids
                            orig_aa = translate(DNAString(codon))
                            new_aa = translate(DNAString(mutated_codon))

                            # Save the amino acid change information
                            aa_changes = rbind(aa_changes, data.frame(Sample_name = sample_name,
                                                                      Chromosome = chrom,
                                                                      Position = snp_pos,
                                                                      Reference_base = ref_base,
                                                                      Alternative_base = alt_base,
                                                                      Exon_start = exon_start,
                                                                      Exon_end = exon_end,
                                                                      Quality = quality,
                                                                      Ref_codon = codon,
                                                                      Alt_codon = mutated_codon,
                                                                      Ref_AA = orig_aa,
                                                                      Alt_AA = new_aa,
                                                                      Total_Depth = total_coverage,
                                                                      Ref_Depth = ref_coverage,
                                                                      Alt_Depth = alt_coverage))
                        }
                    } else {
                        # No overlapping exon found
                        # Save the amino acid change information
                        aa_changes = rbind(aa_changes, data.frame(Sample_name = sample_name,
                                                                  Chromosome = chrom,
                                                                  Position = snp_pos,
                                                                  Reference_base = ref_base,
                                                                  Alternative_base = alt_base,
                                                                  Exon_start = NA,
                                                                  Exon_end = NA,
                                                                  Quality = quality,
                                                                  Ref_codon = NA,
                                                                  Alt_codon = NA,
                                                                  Ref_AA = NA,
                                                                  Alt_AA = NA,
                                                                  Total_Depth = total_coverage,
                                                                  Ref_Depth = ref_coverage,
                                                                  Alt_Depth = alt_coverage))
                    }
                }
            }
        }
    }
}

aa_changes$Frequency = aa_changes$Alt_Depth/aa_changes$Total_Depth
# aa_changes$Sample_name = str_replace_all(aa_changes$Sample_name, ".sorted.bam_SNPs.vcf", "")

write.csv(aa_changes, "/scicore/home/mynaco16/GROUP/analysis_AmpSeq_Daniela_Monica/analysis_Paragon_run1/SNPs_aa_info_v12dec.csv")
