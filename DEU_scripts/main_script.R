# DEXSeq_analysis ####
default_wd <- "/home/dhthutrang/Krebs/Reddy/data"
default_countfolder <- paste(default_wd, "mRNA_seq", "count_files", sep = "/")
default_dexseqfolder <- paste(
  default_wd, "mRNA_seq", "DEXSeq_output",
  sep = "/"
)
default_epigenome1 <- "MCF7_DMSO"
default_epigenome2 <- "MCF7_50nM"
default_referencegenome <- "/home/dhthutrang/ENCODE/refgen/reference_genome.2021_.corrected.gtf"
source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/DEXSeq_analysis.R")

# DEXSeq_analysis 2 ####

source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/utils.R")
