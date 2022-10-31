library(data.table)
library(dplyr)
library(boot)
library(stats)
library(parallel)
library(tidyverse)
library(optparse)
library(doMC)
library(stringr)
doMC::registerDoMC(cores = 17)


option_list <- list(
  make_option(c("--histone_types"),
    type = "character", default = NULL,
    help = "list of histone types separated by comma ','", metavar = "character"
  ), # nolint
  make_option(c("--DEU_file"),
    type = "character", default = NULL,
    help = "Path to DEU output files with format normedcount.csv", metavar = "character"
  ), # nolint
  make_option(c("--DHM_file"),
    type = "character", default = NULL,
    help = "Path to DHM output files with format histone_controlid_treatmentid.bed", metavar = "character"
  ),
  make_option(c("--output_path"),
    type = "character", default = NULL, # nolint
    help = "output path", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)
DEU_file <- "/home/dhthutrang/Krebs/Reddy/data/mRNA_seq/dexseq_output_annotated"
DHM_file <- opt$DHM_file
output_path <- opt$output_path
histone_type_list <- strsplit(opt$histone_types, "_", fixed = T)[[1]]
print(opt)
# ===== PREPARE DATA =====
setwd(output_path)
get_colname <- function(filename_list, option = "his") {
  name <- sapply(filename_list, function(x) strsplit(x, split = "/"))
  name <- sapply(name, function(x) x[length(x)][[1]])
  if (option == "his") {
    name <- sapply(name, function(x) strsplit(x, split = "\\.")[[1]][1])
    name <- sapply(name, function(x) {
      y=strsplit(x, split = "_")[[1]]
      y=paste(y[2:length(y)], collapse = "_")
    })
  } else if (option == "exp") {
    name <- sapply(name, function(x) strsplit(x, split = ".res.csv")[[1]][1])
  }
  # name <- sapply(name, function(x) paste(sort(strsplit(x, split = "_")[[1]]), collapse = "_"))
  names(name) <- NULL
  print(name)
  return(name)
}

# ===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====
print("===== PREPARE EXP FILE (1 FOR ALL HIS TYPES) =====")
all_pairs.exp <- list.files(DEU_file, full.names = TRUE, pattern = "grouped.bed")
print(all_pairs.exp)
get_all_pairs.exp <- function(all_pairs.exp) {
  pair.exp <- fread(all_pairs.exp[[1]])
  id <- as.data.frame(do.call(rbind, lapply(pair.exp$V9, function(x) strsplit(x, split = '"', fixed = T)[[1]][c(2, 6)])))
  colnames(id) <- c("gene", "exon")
  pair.exp <- pair.exp %>%
    dplyr::mutate(
      gene_id = id$gene, exon_id = id$exon,
      MCF7_DMSO_MCF7_50nM = V10
    ) %>%
    select(gene_id, exon_id, MCF7_DMSO_MCF7_50nM) %>%
    dplyr::na_if(., -Inf) %>%
    dplyr::na_if(., Inf)
    
  print(dim(pair.exp))
  print(head(pair.exp, 10))
  # saveRDS(pair.his_list, paste('pair.his_list_', his, '.RDS', sep=''))
  return(pair.exp)
}

all_pairs.exp <- get_all_pairs.exp(all_pairs.exp)
saveRDS(all_pairs.exp, "all_pairs.exp.RDS")
all_pairs.exp <- readRDS("all_pairs.exp.RDS")
head(all_pairs.exp)

# # ===== PREPARE HIS FILE (6 TOTAL) =====
print("===== PREPARE HIS FILE (6 TOTAL) =====")
# his_id = read.csv("flank_id.2021.txt", sep='\t', header = FALSE)
get_all_pairs.his <- function(all_pairs.his, his) {
  pair.his_list <- vector("list", length(all_pairs.his))
  for (i in 1:length(all_pairs.his)) {
    # for (i in 1:1){
    print(paste("Pair: ", i, sep = ""))
    pair.his <- all_pairs.his[[i]]
    pair.his <- fread(pair.his)
    id <- as.data.frame(do.call(rbind, lapply(pair.his$V9, function(x) strsplit(x, split = '"', fixed = T)[[1]][c(2, 6)])))
    colnames(id) <- c("gene", "exon")
    pair.his <- pair.his %>%
      dplyr::mutate(
        gene = id$gene, exon = id$exon, type = V3,
        p_val = as.numeric(as.character(V11)),
        m_val = dplyr::if_else(p_val <= 0.05,
          true = abs(as.numeric(as.character(V10))), false = 0
        )
      ) %>%
      dplyr::select(gene, exon, m_val) %>%
      dplyr::na_if(., -Inf) %>%
      dplyr::na_if(., Inf)
    if (i == 1) pair.his_id <- pair.his[, c("gene", "exon")]
    pair.his <- pair.his %>% dplyr::select(-gene, -exon)
    pair.his_list[[i]] <- pair.his
  }
  lapply(pair.his_list, function(x) print(dim(x)))
  pair.his_list <- as.data.frame(cbind(pair.his_id, as.data.frame(do.call(cbind, pair.his_list))))
  print(dim(pair.his_list))
  print(head(pair.his_list))
  # saveRDS(pair.his_list, paste('pair.his_list_', his, '.RDS', sep=''))
  return(pair.his_list)
}

get_all_pairs.his_list <- function(histone_type_list) {
  all_pairs.his_list <- vector("list", length(histone_type_list))
  for (j in 1:length(histone_type_list)) {
    his <- histone_type_list[j]
    print(his)
    all_pairs.his <- list.files(DHM_file, pattern = ".bed", full.names = TRUE)
    all_pairs.his <- all_pairs.his[grep(his, all_pairs.his)]
    print(all_pairs.his)
    colname_his <- c("gene_id", "exon_id", get_colname(all_pairs.his, "his"))
    all_pairs.his.sig <- get_all_pairs.his(all_pairs.his, his)
    colnames(all_pairs.his.sig) <- colname_his

    all_pairs.his_list[[j]] <- as.data.table(all_pairs.his.sig)
  }
  return(all_pairs.his_list)
}
filter_all_his_list <- function(his_list, histone_type_list, filter_genes_path) {
  all_filtered_df <- vector("list")
  for (i in 1:length(histone_type_list)) {
    histone <- histone_type_list[i]
    his_df <- his_list[[i]]
    print(i)
    print(histone)
    all_filtered_df[[i]] <- filter_genes(df = his_df, filter_genes_path = filter_genes_path, filter = histone)
  }
  names(all_filtered_df) <- histone_type_list
  return(all_filtered_df)
}

all_pairs.his_list = get_all_pairs.his_list(histone_type_list)
saveRDS(all_pairs.his_list, "all_pairs.his_list.RDS")
all_pairs.his_list <- readRDS("all_pairs.his_list.RDS")
# names(all_pairs.his_list) <- histone_type_list


# ===== CORRELATION WITH RANDOMIZATION =====
# ------------ Execute analysis ------------
# included_genes = fread("gene_id.txt", header = FALSE)
# included_genes = unique(all_pairs.exp_flt$gene_id)
p_value_calculator <- function(r, nrow) {
  P <- r * sqrt(nrow - 2) / sqrt(1 - r * r)
  P <- 2 * pt(-abs(P), nrow - 2)
  return(P)
}
pearcor_p <- function(exp, his) {
  if (length(unique(exp)) > 1 & length(unique(his)) > 1) {
    p_val <- p_value_calculator(cor(exp, his, method = "pearson"), nrow = length(exp))
    # p_val = p.adjust(p_val, method = "fdr", n=ncol(included_genes))
    return(p_val)
  } else {
    return(NA)
  }
}
pearcor_r <- function(exp, his, n_points) {
  df <- as.data.frame(cbind(exp, his))
  n_sep_point <- nrow(unique(df))
  if (0 %in% apply(df, 1, unique)) n_sep_point <- n_sep_point - 1
  if (n_sep_point >= 3 & length(unique(exp)) > 1 & length(unique(his)) >= n_points) {
    r_val <- cor(exp, his, method = "pearson")
    return(r_val)
  } else {
    return(NA)
  }
}
analyze_array <- function(all_pairs.exp, all_pairs.his, option = "p", n_points, included_genes) {
  all_res_pair <- vector("list", ncol(all_pairs.exp) - 2)
  subset_name <- colnames(all_pairs.his)
  print(subset_name)
  all_pairs.exp_subset <- all_pairs.exp[, ..subset_name]

  # for (i in 1:n_pairs){
  all_res_pair <- foreach(i = 1:(length(subset_name) - 2), .combine = "c", .packages = c("dplyr")) %dopar% {
    print(paste("Pair: ", i, sep = ""))
    exp <- all_pairs.exp_subset[[i + 2]]
    his <- all_pairs.his[[i + 2]]
    data_table <- as.data.table(cbind(exp, his))
    head(data_table)

    if (option == "p") {
      res_table <- data_table %>%
        dplyr::group_by(all_pairs.exp$gene_id) %>%
        dplyr::summarise(res = pearcor_p(exp, his)) %>%
        dplyr::select(res)
    } else if (option == "r") {
      res_table <- data_table %>%
        dplyr::group_by(all_pairs.exp$gene_id) %>%
        dplyr::summarise(res = pearcor_r(exp, his, n_points)) %>%
        dplyr::select(res)
    }
    # all_res_pair[[i]] = res_table
  }
  all_res_pair <- as.data.table(all_res_pair)
  all_res_pair <- cbind(included_genes, all_res_pair)
  colnames(all_res_pair) <- c("gene_id", subset_name[3:length(subset_name)])
  print(head(all_res_pair))
  return(as.data.frame(all_res_pair))
}
analyze_array_list <- function(all_pairs.exp, all_pairs.his_list, method = "p", n_points = 2) {
  all_res_list <- vector("list", length(histone_type_list) - 1)
  included_genes <- unique(all_pairs.exp$gene_id)
  for (j in 1:length(histone_type_list)) {
    print(paste("Histone: ", histone_type_list[j], sep = ""))
    all_pairs.his <- all_pairs.his_list[[j]]
    if (method == "p") {
      all_res_pair <- analyze_array(all_pairs.exp, all_pairs.his, option = "p", included_genes = included_genes)
    } else if (method == "r") {
      all_res_pair <- analyze_array(all_pairs.exp, all_pairs.his, option = "r", n_points = n_points, included_genes = included_genes)
    }
    all_res_list[[j]] <- all_res_pair
  }
  return(all_res_list)
}

print("Pearsons-p correlation")
all_res_list.pearcor_p <- analyze_array_list(all_pairs.exp, all_pairs.his_list, method = "p")
names(all_res_list.pearcor_p) = histone_type_list
saveRDS(all_res_list.pearcor_p, "all_res_list.pearcor_p.RDS")

all_res_list.pearcor_r <- analyze_array_list(all_pairs.exp, all_pairs.his_list, method = "r")
names(all_res_list.pearcor_r) = histone_type_list
saveRDS(all_res_list.pearcor_r, "all_res_list.pearcor_r.RDS")

head(r_val[[1]][order(r_val[[1]][2]),])

r_val = readRDS("/home/dhthutrang/Krebs/Reddy/data/cor_analysis/all_res_list.pearcor_r.RDS")
r_val = lapply(r_val, function(x) x[!is.na(x[[2]]) & abs(x[[2]]) > 0.3,])
r_val = Reduce(union, r_val)
r_val
p_val = readRDS("/home/dhthutrang/Krebs/Reddy/data/cor_analysis/all_res_list.pearcor_p.RDS")
p_val = lapply(p_val, function(x) x[!is.na(x[[2]]) & x[[2]] <= 0.05, ]) # nolint
p_val = Reduce(union, p_val)
p_val
temp = intersect(p_val, r_val)
H3K27ac_res = r_val[['H3K27ac']][r_val[['H3K27ac']]$gene_id %in% intersect(r_val[['H3K27ac']]$gene_id, p_val[['H3K27ac']]$gene_id), ]
H3K9ac_res = r_val[['H3K9ac']][r_val[['H3K9ac']]$gene_id %in% intersect(r_val[['H3K9ac']]$gene_id, p_val[['H3K9ac']]$gene_id), ]
# lapply(p_val, head)
p_val[[2]][p_val[[2]][1] == 'C1QBP']

temp = unlist(sapply(temp, function(x) strsplit(x, "+", fixed = TRUE)[[1]]))
paste(temp[order(temp)], collapse = ", ")
print(H3K27ac_res[order(H3K27ac_res[2]),])
print(H3K9ac_res[order(H3K9ac_res[2]),])
