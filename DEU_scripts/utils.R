library(data.table, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(reshape2, quietly = TRUE)

library("optparse", quietly = TRUE)
library("BiocParallel")


original_wd <- getwd()

option_list <- list(
    make_option(c("-i", "--countfolder"),
        type = "character", default = default_countfolder,
        help = "path to folder of counts", metavar = "character"
    ),
    make_option(c("-o", "--default_dexseqfolder"),
        type = "character", default = default_dexseqfolder,
        help = "path to folder of result", metavar = "character"
    ),
    make_option(c("-a", "--epigenome1"),
        type = "character", default = default_epigenome1,
        help = "ID of first epigenome", metavar = "character"
    ),
    make_option(c("-b", "--epigenome2"),
        type = "character", default = default_epigenome2,
        help = "ID of second epigenome", metavar = "character"
    ),
    make_option(c("-g", "--referencegenome"),
        type = "character", default = default_referencegenome,
        help = "path to flattened reference genome", metavar = "character"
    ),
    make_option(c("-n", "--numcores"),
        type = "integer", default = 1,
        help = "number of processing cores", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
dexseqfolder <- opt$default_dexseqfolder
dir.create(file.path(dexseqfolder, "analysis"), showWarnings = FALSE)
dexseq_analysis_folder <- paste(dexseqfolder, "analysis", sep = "/")

get_deu_genes <- function(
                          fold_change_index = 7,
                          fold_change_threshold = 2.0,
                          adj_p_index = 6,
                          adj_p_threshold = 0.0001) {
    print("----> Getting DEU genes")
    setwd(dexseqfolder)
    dir_files <- list.files(pattern = "*.csv")
    deu_genes_list <- vector("list", length(dir_files))
    for (i in 1:length(dir_files)) {
        print(dir_files[i])
        data <- data.table::data.table(read.csv(dir_files[i], sep = "\t"))
        sig_exon <- data[
            abs(data[[fold_change_index]]) >= fold_change_threshold &
                data[[adj_p_index]] <= 0.0001,
        ]
        sig_exon <- sig_exon %>%
            dplyr::mutate(gene_id = paste(groupID, featureID, sep = ";"))
        deu_genes_list[[i]] <- sig_exon$gene_id
    }
    deu_genes_list <- Reduce(c, deu_genes_list)
    saveRDS(
        deu_genes_list,
        paste("analysis", "deu_genes_list.RDS", sep = "/")
    )
    setwd(original_wd)
    print("----> Finished")
    return(deu_genes_list)
}

summarize_deu_genes_info <- function(deu_genes_list) {
    exon_freq <- table(deu_genes_list)
    exon_freq <- exon_freq[order(exon_freq, decreasing = TRUE)]
    gene_freq <- table(
        unlist(lapply(
            deu_genes_list,
            function(x) strsplit(x, ";", fixed = TRUE)[[1]][1]
        ))
    )
    gene_freq <- gene_freq[order(gene_freq, decreasing = TRUE)]
    paste(names(gene_freq)[gene_freq > 3], collapse = ", ")
}

deu_genes_list <- get_deu_genes()
summarize_deu_genes_info()