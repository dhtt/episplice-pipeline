start_time <- Sys.time()
# ===== LOAD PACKAGES ======
library("stringr", quietly = TRUE)
library("data.table", quietly = TRUE)
library("dplyr", quietly = TRUE)
library("DEXSeq", quietly = TRUE)
library("optparse", quietly = TRUE)
library("BiocParallel")

option_list <- list(
  make_option(c("-i", "--countfolder"),
    type = "character", default = NULL,
    help = "path to folder of counts", metavar = "character"
  ),
  make_option(c("-o", "--resultfolder"),
    type = "character", default = NULL,
    help = "path to folder of result", metavar = "character"
  ),
  make_option(c("-a", "--epigenome1"),
    type = "character", default = NULL,
    help = "ID of first epigenome", metavar = "character"
  ),
  make_option(c("-b", "--epigenome2"),
    type = "character", default = NULL,
    help = "ID of second epigenome", metavar = "character"
  ),
  make_option(c("-g", "--referencegenome"),
    type = "character", default = NULL,
    help = "path to flattened reference genome", metavar = "character"
  ),
  make_option(c("-n", "--numcores"),
    type = "integer", default = 1,
    help = "number of processing cores", metavar = "character"
  )
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ===== PREPARE DATA =====
setwd(opt$countfolder)
in_dir <- normalizePath(getwd())
epi_id1 <- opt$epigenome1
epi_id2 <- opt$epigenome2
pair <- paste(
  paste("^", epi_id1, ".*.txt$", sep = ""),
  paste("^", epi_id2, ".*.txt$", sep = ""), sep = "|"
  )
count_files <- list.files(in_dir, pattern = pair, full.names = TRUE)
file_names <- as.data.table(str_split_fixed(basename(count_files), "\\_", 3))
gtf_files <- opt$referencegenome
cores <- MulticoreParam(opt$numcores)

print(paste("Working folder: ", opt$countfolder, sep = ""))
print(paste("Epigenomes: ", epi_id1, epi_id2, sep = ""))
print(paste("Count files: ", count_files, sep = ""))
print(paste("GTF files: ", gtf_files, sep = ""))

print(file_names)
sample_table <- data.frame(
  #   row.names = c(paste(c(file_names$V1, file_names$V2, file_names$V3), sep="_")),
  condition = c(file_names$V2)
)
print("Sample table: ")
print(sample_table)

# ===== RUN DEXSEQ =====
print("---> Inputting to DEXSeq")
dxd <- DEXSeqDataSetFromHTSeq(
  count_files,
  sampleData = sample_table,
  design = ~ sample + exon + condition:exon,
  flattenedfile = normalizePath(gtf_files)
)

print("Getting DEXSeq result")
dxd_res <- DEXSeq(dxd, quiet = FALSE, BPPARAM = cores)

# ===== SAVING RESULTS =====
print("Saving DEXSeq normalized counts")
dxd_count <- data.frame(cbind(
  dxd_res[c(1, 2)],
  counts(dxd_res, normalized = TRUE)
  ))
colnames(dxd_count) <- c(
  "groupID",
  "featureID",
  paste(file_names$V1, file_names$V2, sep = "_")
)
normedcount_path <- paste(
  opt$resultfolder,
  paste(
    paste(
      epi_id1, epi_id2, sep = "_"),
      "normedcount.csv", sep = "_"),
      sep = "/"
      )
write.table(
  dxd_count, normedcount_path,
  quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE
)

print("Saving DEXSeq result")
r_data_path <- paste(
  opt$resultfolder,
  paste(
    paste(epi_id1, epi_id2, sep = "_"),
    "RData", sep = "."),
    sep = "/"
    )
save(dxd_res, file = r_data_path)

result_path <- paste(
  opt$resultfolder,
  paste(
    paste(epi_id1, epi_id2, sep = "_"),
    "res.csv", sep = "_"),
    sep = "/"
    )
write.table(as.data.frame(dxd_res[c(1, 2, 3, 5, 6, 7, 10)]), result_path,
  quote = FALSE, sep = "\t", dec = ".", row.names = FALSE, col.names = TRUE
)

print("---> Exporting HTML DEXSeq result")
html_path <- paste(
  opt$resultfolder,
  paste(paste(epi_id1, epi_id2, sep = "_"), "html", sep = "_"),
  sep = "/"
  )
DEXSeqHTML(dxd_res,
  path = html_path,
  FDR = 0.05, color = c("#FF000080", "#0000FF80"),
  BPPARAM = cores
)
print("===> FINISHED!")

end_time <- Sys.time()
print(paste("Total time:", end_time - start_time, sep = " "))
