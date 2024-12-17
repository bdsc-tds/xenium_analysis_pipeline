log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Annotate Xenium from Chromium using XGBoost

library(Seurat)
library(dplyr)
library(xgboost)

# Load common reference-based parameters
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_header.R")
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_xgboost_classifications.R")

# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'counts' 
test_layer <- 'counts' 

# Load reference and query data
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_load_data.R") # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

CELL_MIN_INSTANCE <- snakemake@params[["cell_min_instance"]]

# Generate reference object
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_generate_reference_obj.R")

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

# Predict probabilities
xe_pred_prob <- predict(regressions_classification$fit, GetAssayData(xe, assay = xe_assay, layer = test_layer)[xe_chrom_common_genes,] %>% Matrix::t())
xe_pred_prob <- matrix(xe_pred_prob, ncol = length(Y.ref_map), byrow = TRUE)
colnames(xe_pred_prob) <- Y.ref_map %>% as.vector()
rownames(xe_pred_prob) <- colnames(xe)

# Convert probabilities to class labels
labels <- max.col(xe_pred_prob) - 1  # Subtract 1 because `max.col` returns 1-based index
labels <- Y.ref_map[as.character(labels)] # convert numbers to cell types 
names(labels) <- colnames(xe)

# Save annotation
saveRDS(regressions_classification, snakemake@output[[1]])
write.csv(labels, snakemake@output[[2]])
write.csv(xe_pred_prob, snakemake@output[[3]]) 


