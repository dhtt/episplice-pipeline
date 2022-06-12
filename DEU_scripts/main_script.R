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
# source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/DEXSeq_analysis.R")

# DEXSeq analysis report ####
# source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/utils.R")
# deu_genes_list <- get_deu_genes()
# summarize_deu_genes_info()

# DEXSeq analysis report ####
default_epigenome1 <- NULL
default_epigenome2 <- NULL
default_gene_id <- NULL
source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/utils.R", encoding = "UTF-8")
get_dexseq_plot(gene_id = "STIM1")
get_dexseq_plot(epi_id1 = "MCF7_DMSO", epi_id2 = "MCF7_50nM", gene_id = "STIM1")
