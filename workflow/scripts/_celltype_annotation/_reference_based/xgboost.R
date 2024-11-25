# Annotate Xenium from Chromium using XGBoost

library(Seurat)
library(dplyr)
library(xgboost)

# Load common reference-based parameters
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_header.R")) 


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'counts' 
test_layer <- 'counts' 

# Load reference and query data
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_load_data.R")) # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

# Generate reference object
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_generate_reference_obj.R")) 

ref_labels <- chrom@meta.data %>% pull(annotation_level) %>% as.vector()

xe_chrom_common_genes <- intersect(rownames(xe), rownames(chrom))

# Train xgboost
# Get xgboost params
nrounds <- snakemake@params[["nrounds"]] %||% 1000 # can be xgboost-specific param
eta <- snakemake@params[["eta"]] %||% 0.3 # can be xgboost-specific param

# Transform labels to integers for compatibility with xgboost
ref_labels_factor <- factor(chrom@meta.data %>% pull(annotation_level))
Y.ref <- as.integer(ref_labels_factor) - 1
Y.ref_map <- levels(ref_labels_factor)
names(Y.ref_map) <- seq_along(Y.ref_map) - 1

# Train
regressions_classification <- run_regressions_classification_fixed_features(
  X = GetAssayData(chrom, assay = ref_assay, layer = ref_layer)[xe_chrom_common_genes,] %>% Matrix::t(),
  Y = Y.ref, 
  do.cross.val = FALSE,
  nrounds = nrounds, 
  eta = eta, 
  nthread = parallel::detectCores(),
  verbose = FALSE
)

# Annotate with XGBoost 
labels   <- predict(regressions_classification$fit, GetAssayData(xe, assay = xe_assay, layer = test_layer)[xe_chrom_common_genes,] %>% Matrix::t())
labels   <- Y.ref_map[as.character(labels)]
names(labels) <- colnames(xe)

# Save annotation
saveRDS(regressions_classification, snakemake@output[[1]])
write.csv(labels, snakemake@output[[2]])
#write.csv(scores, snakemake@output[[3]]) # no score-like values in xgboost? 


