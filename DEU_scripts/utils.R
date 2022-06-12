library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("ggplot2", quietly = TRUE)
library("reshape2", quietly = TRUE)
library("DEXSeq", quietly = TRUE)

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
    make_option(c("-G", "--referencegenome"),
        type = "character", default = default_referencegenome,
        help = "path to flattened reference genome", metavar = "character"
    ),
    make_option(c("-n", "--numcores"),
        type = "integer", default = 1,
        help = "number of processing cores", metavar = "character"
    ),
    make_option(c("-g", "--gene_id"),
        type = "character", default = NULL,
        help = "chosen gene for plotting", metavar = "character"
    )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
dexseqfolder <- opt$default_dexseqfolder
dir.create(file.path(dexseqfolder, "analysis"), showWarnings = FALSE)
dexseq_analysis_folder <- paste(dexseqfolder, "analysis", sep = "/")
dexseq_normedcount_path <- paste(dexseqfolder, "csv", sep = "/")
dexseq_result_path <- paste(dexseqfolder, "csv", sep = "/")
dexseq_r_data_path <- paste(dexseqfolder, "r_data", sep = "/")
dexseq_html_path <- paste(dexseqfolder, "html", sep = "/")
dexseq_plot_path <- paste(dexseqfolder, "plot", sep = "/")

epi_id1 <- opt$epigenome1
epi_id2 <- opt$epigenome2
gene_id <- opt$gene_id
image_height <- opt$height
image_width <- opt$width
image_res <- opt$resolution

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

get_dexseq_plot <- function(
                            epi_id1 = epi_id1, epi_id2 = epi_id2,
                            gene_id = gene_id, image_res = 400,
                            image_height = 8, image_width = 5) {
    print("----> Plotting")
    setwd(dexseq_r_data_path)
    r_data_filename <- paste(
        paste(paste(epi_id1, epi_id2, sep = "_"), "*.RData", sep = "."),
        paste(paste(epi_id2, epi_id1, sep = "_"), "*.RData", sep = "."),
        sep = "|"
    )
    r_data_files <- list.files(
        normalizePath(getwd()),
        pattern = r_data_filename, full.names = FALSE
    )[1]
    print(r_data_files)

    load(r_data_files)
    plot_name <- paste(
        paste(dexseq_plot_path, r_data_files, sep = "/"),
        gene_id, "tiff",
        sep = "."
    )
    tiff(
        plot_name,
        width = image_width, height = image_height,
        units = "in", res = image_res
    )
    plotDEXSeq(
        dxd.res, gene_id,
        legend = TRUE, FDR = 0.05,
        norCounts = TRUE, splicing = TRUE, expression = TRUE,
        color = c("#EE442F", "#63ACBE"), lwd = 1.8
    )
    dev.off()
    setwd(original_wd)
    print("----> Finished")
}
