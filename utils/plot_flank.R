library(Gviz)
library(GenomicRanges)
library(rtracklayer)
library(stringr)
library(data.table)
library(dplyr)
library(optparse)
library(ggplot2)

# option_list <- list(
#   make_option(c("--example_path"),
#     type = "character", default = NULL,
#     help = "Path to example gene','", metavar = "character"
#   ) # nolint
# #   make_option(c("--DEU_file"),
# #     type = "character", default = NULL,
# #     help = "Path to DEU output files with format normedcount.csv", metavar = "character"
# #   ), # nolint
# #   make_option(c("--DHM_file"),
# #     type = "character", default = NULL,
# #     help = "Path to DHM output files with format histone_controlid_treatmentid.bed", metavar = "character"
# #   ),
#   make_option(c("--output_path"),
#     type = "character", default = NULL, # nolint
#     help = "output path", metavar = "character"
#   )
# )

# opt_parser <- OptionParser(option_list = option_list)
# opt <- parse_args(opt_parser)
# example_path <- opt$example_path
# output_path <- opt$output_path

make_equidist <- function(gff_df){
    start <- min(gff_df$start)
    end <- max(gff_df$end)

    new_coords <- seq(start, end, abs(start - end)/(nrow(gff_df)*2-1))
    new_start <- new_coords[seq(1,length(new_coords), 2)]
    new_end <- new_coords[seq(2,length(new_coords), 2)]
    # print(new_coords)
    # print(DataFrame(cbind(new_start, new_end)))

    gff_df_equi <- gff_df
    gff_df_equi$start_equi = new_start
    gff_df_equi$end_equi = new_end
    return(gff_df_equi)
}

example_path <- "/home/dhthutrang/Krebs/episplice-pipeline/utils/result/CAB39L"
output_path = example_path

col.highlight <- "#fffa70"
fontsize.title <- 17 #normally 17
type.color <- c("#EE442F", "#63ACBE")
margin <- c(20, 20)
background.title <- "#5D2ED2"
col.title <- "white"

all_files <- list.files(normalizePath(example_path), full.names = TRUE)

#============== Plot gene region ==============
collapsed_flanks <- import.gff(all_files[grep("FLANK.txt", all_files)]) %>%
    as("GRanges")
collapsed_flanks_equi <- make_equidist(collapsed_flanks)

gen <- 'hg19'
chr <- as.character(unique(seqnames(collapsed_flanks)))[[1]]
gene_name <- collapsed_flanks$gene_id[1]
reverseStrand <- TRUE

itrack <- IdeogramTrack(genome = gen, chromosome = chr, name = gene_name,
    fontsize.title = fontsize.title
)
gtrack <- GenomeAxisTrack(labelPos = "above", 
fontsize.title = fontsize.title
)

collapsed_flanks_track <- GeneRegionTrack(collapsed_flanks,
                           genome = gen, chromosome = chr,
                           name = "Transcripts", stacking = 'hide',
                           fontsize.title = fontsize.title,
                           col = "black", fill="#85c0f9",
                           cex.group = 0.5
)
collapsed_flanks_track_equi <- GeneRegionTrack(collapsed_flanks_equi,
                           genome = gen, chromosome = chr,
                           name = "Flattened", 
                           fontsize.title = fontsize.title,
                           col = "black", fill="#85c0f9",
                           cex.group = 0.5
)

#============== Plot gene ===========
all_tracks <- list(itrack, gtrack, collapsed_flanks_track, collapsed_flanks_track_equi)
tiff(paste(output_path, paste(gene_name, 'tiff', sep='.'), sep='/'), units="in", width=12, height=5, res=300)
plotTracks(all_tracks,
           extend.left = 5000, extend.right = 500,
           background.title = background.title,
           fontsize.title = fontsize.title,
           # fontsize = fontsize.title,
           col.title = col.title,
           margin = margin, 
           reverseStrand = reverseStrand,
           stackHeight = 0.75, cex = 0.6, lwd=1.5, frame=TRUE, cex.axis= 0.3, 
           lineheight = 0.25, cex.legend= 0.75, fontsize.legend = fontsize.title)
dev.off()
#=================
collapsed_flanks <- import.gff(all_files[grep("FLANK.txt", all_files)]) %>%
    as('data.frame') 
# collapsed_flanks$type = 'collapsed_flanks'
collapsed_flanks$floor = 1.5
collapsed_flanks$ceiling = 2.5
collapsed_flanks_equi <- make_equidist(collapsed_flanks)
# collapsed_flanks_equi$type = 'collapsed_flanks_equi'
collapsed_flanks_equi$floor_equi = 0
collapsed_flanks_equi$ceiling_equi = 1

DEU <- fread(all_files[grep("DEU", all_files)]) %>%
    dplyr::mutate(
        exon_no = V2,
        DEU = if_else(is.na(V4), 0, V4),
        FDR = if_else(is.na(V6), 1, V6)
    ) %>%
    select(exon_no, DEU, FDR)
dim(collapsed_flanks_equi)

data <- cbind(collapsed_flanks_equi)


plot <- ggplot(data=data) +
    geom_rect(aes(xmin=start, xmax=end, ymin=floor, ymax=ceiling), fill='#0095eb') +
    geom_rect(aes(xmin=start_equi, xmax=end_equi, ymin=floor_equi, ymax=ceiling_equi), fill='#7400b8') +
    geom_segment(data=data, aes(xend=(start+end)/2, x=(start_equi + end_equi)/2, 
        y=ceiling_equi, yend=floor), col='grey')


tiff(paste(output_path, paste(gene_name, 'ggplot.tiff', sep='.'), sep='/'), units="in", width=12, height=5, res=300)
plot
dev.off()


exp = readRDS('/home/dhthutrang/Krebs/Reddy/data/cor_analysis/all_pairs.exp.RDS')
exp[exp$gene_id == "CAB39L",]
his = readRDS('/home/dhthutrang/Krebs/Reddy/data/cor_analysis/all_pairs.his_list.RDS')
his = his[[1]]
his[his$gene_id == "CAB39L",]
cor(exp[exp$gene_id == "CAB39L",3], his[his$gene_id == "CAB39L",3])
length(his[his$gene_id == "CAB39L",3])
length(exp[exp$gene_id == "CAB39L",3])
