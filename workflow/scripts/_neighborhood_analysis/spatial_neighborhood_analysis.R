library(Seurat)
library(spacexr)
library(arrow)
library(dplyr)
library(puRCTD)

xe <- readRDS(snakemake@input[["xe"]])
rctd <- readRDS(snakemake@input[["post_processed_rctd"]])

### Compute Spatial network & conduct analysis ###
sp_nw <- build_spatial_network(
  xe, 
  dims = 1:2, 
  DO_prune = snakemake@params[["DO_prune"]], 
  rad_pruning = snakemake@params[["pruning_radius"]],
  k_knn = snakemake@params[["k_knn"]]
  )

sp_nw <- add_spatial_metric(spatial_neighborhood = sp_nw, rctd = rctd)
sp_neigh_df <- neighborhood_analysis_to_metadata(sp_nw)

# Output spatial scores
write_parquet(sp_neigh_df, snakemake@output[["spatial_neighborhood_scores"]])
