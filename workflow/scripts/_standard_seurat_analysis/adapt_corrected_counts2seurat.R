log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(Seurat)

xe <- readRDS(snakemake@input[["raw_obj"]])
corrected_counts <- Read10X_h5(snakemake@input[["corrected_counts"]])

# Substitute original reads with corrected ones @Mariia

# Save object
saveRDS(
  xe, 
  file = file.path(snakemake@output[[1]])
)
