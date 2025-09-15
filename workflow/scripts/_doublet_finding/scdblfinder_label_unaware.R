log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)
library(spacexr)
library(SingleCellExperiment)
library(scater)
library(scDblFinder)
library(arrow)
library(dplyr)


xe <- readRDS(snakemake@input[["xe"]])
rctd <- readRDS(snakemake@input[["post_processed_rctd"]])

################################## PURIFICATION ############################################
message("Adding metadata to Xenium obj... \n")
rctd_df <- rctd@results$results_df
xe <- AddMetaData(xe, rctd_df)
DefaultAssay(xe) <- "Xenium"

message("Converting Seurat to SCE ... \n")
xe_sce <- as.SingleCellExperiment(xe)
message("Done converting Seurat to SCE\n")

message("Running label-UNaware scDblFinder...\n")
set.seed(123)
Nclusters <- xe_sce$first_type %>% unique() %>% length()
xe_sce <- scDblFinder(xe_sce, clusters = Nclusters, dbr.sd = 1) # using large sd as suggested by vignette when doublet rate is unknown
message("Done running label-UNaware scDblFinder...\n")

res <- xe_sce@colData

message("Saving scDblFinder output ...\n")
# Output purified counts
write_parquet(as.data.frame(res), snakemake@output[["doublet_scores"]])
message("Done saving scDblFinder output \n")