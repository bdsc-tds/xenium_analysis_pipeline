log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)
library(spacexr)
library(arrow)
library(dplyr)
library(SPLIT)

options(future.globals.maxSize = snakemake@params[["future_globals_maxSize"]])

xe <- readRDS(snakemake@input[["xe"]])
rctd <- readRDS(snakemake@input[["post_processed_rctd"]])

### Compute Transcriptomic network & conduct analysis ###
tr_nw <- build_transcriptomics_network(
  xe,
  DO_prune = snakemake@params[["DO_prune"]],
  k_knn = snakemake@params[["k_knn"]]
)
tr_nw <- add_transcriptomics_metric(transcriptomics_neighborhood = tr_nw, rctd = rctd) 
tr_neigh_df <- neighborhood_analysis_to_metadata(tr_nw)

# Output transcriptomic scores
write_parquet(tr_neigh_df, snakemake@output[["transcriptomic_neighborhood_scores"]])
