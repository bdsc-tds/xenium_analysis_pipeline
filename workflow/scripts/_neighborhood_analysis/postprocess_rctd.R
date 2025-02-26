library(spacexr)
library(arrow)
library(dplyr)
library(puRCTD)

rctd <- readRDS(snakemake@input[["rctd_result"]])

### Post RCTD -- correct rctd scores and compute additional
rctd <- puRCTD::run_post_process_RCTD(rctd = rctd)

### Output post-processed RCTD
results_df <- rctd@results$results_df
results_df$cell_id <- rownames(results_df)
results_df <- results_df %>% select(cell_id, everything())

write_parquet(results_df, snakemake@output[["post_processed_rctd_df"]])

saveRDS(rctd, file = snakemake@output[["post_processed_rctd"]])
