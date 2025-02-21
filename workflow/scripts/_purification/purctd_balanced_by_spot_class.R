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
purified_counts <- Read10X_h5(snakemake@input[["fully_purified_counts"]])
purified_counts_metadata <- read_parquet(snakemake@input[["fully_purified_counts_metadata"]])

###### Balancing purification level ######

# Create purified Seurat
xe_purified <- CreateSeuratObject(counts = purified_counts, meta.data = purified_counts_metadata)
rm("purified_counts")

### Run `spot_class_based` purification (ie, purify doublets and do not touch singlets)
### Done by balancing raw and purified datasets
xe <- AddMetaData(xe, rctd@results$results_df)
xe_balanced_spot_class <- balance_raw_and_purified_data_by_spot_class(
  xe_raw = xe,
  xe_purified = xe_purified,
  spot_class_key = "spot_class"
)
# Output spot-class-based purified counts
write10xCounts(
  path = snakemake@output[["spot_class_purified_counts"]], 
  x = GetAssayData(xe_balanced_spot_class, assay = "RNA", layer = "counts")
)
write_parquet(
  xe_balanced_spot_class@meta.data %>% select(all_of(colnames(purified_counts_metadata))), 
  snakemake@output[["spot_class_purified_counts_metadata"]]
)
