# DEXSeq_analysis ####
default_wd <- "/home/dhthutrang/Krebs/Shen"
default_countfolder <- paste(default_wd, "mRNA_seq", "count_files", sep = "/")
default_dexseqfolder <- paste(
  default_wd, "mRNA_seq", "DEXSeq_output",
  sep = "/"
)
default_epigenome1 <- "MCF7_parental"
default_epigenome2 <- "MCF7_KO"
default_referencegenome <- "/home/dhthutrang/ENCODE/refgen/reference_genome.2021_.corrected.gtf" # nolint
# source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/DEXSeq_analysis.R") # nolint

# DEXSeq analysis report ####
source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/utils.R", encoding = "UTF-8") # nolint
# deu_genes_list <- get_deu_genes(adj_p_threshold = 0.01)
# summarize_deu_genes_info(deu_genes_list)

# DEXSeq analysis report ####
default_epigenome1 <- NULL
default_epigenome2 <- NULL
default_gene_id <- NULL
source("/home/dhthutrang/Krebs/episplice-pipeline/DEU_scripts/utils.R", encoding = "UTF-8") # nolint
# get_dexseq_plot(epi_id1 = "MCF7_parental", epi_id2 = "MCF7_KO", gene_id = "STIM1") # nolint

library(RCurl)
x = getURL("https://raw.githubusercontent.com/nf-core/test-datasets/chipseq/design_full.csv")