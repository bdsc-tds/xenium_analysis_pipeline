library(Seurat)
library(spacexr)
library(arrow)
library(dplyr)

if (!requireNamespace("puRCTD", quietly = TRUE)){
  remotes::install_git("git@github.com:bdsc-tds/puRCTD.git")
}
library(puRCTD)

snakemake@source("../../../scripts/utils/readwrite.R") # @Senbai, this one is from `jon_count_correction` branch


xe <- readRDS(snakemake@input[["xe"]])
rctd <- readRDS(snakemake@input[["post_processed_rctd"]])

################################## PURIFICATION ############################################
### Run full purification (ie., purify all cells except for highly confident singlets (ie., those that do not have the second cell type)
res_purification <- puRCTD::purify(
  counts = GetAssayData(xe, assay = 'Xenium', layer = 'counts'),
  rctd = rctd,
  DO_purify_singlets = TRUE
)
# Output purified counts
write10xCounts(path = snakemake@output[["fully_purified_counts"]], x = res_purification$purified_counts)
write_parquet(res_purification$cell_meta, snakemake@output[["fully_purified_counts_metadata"]])

