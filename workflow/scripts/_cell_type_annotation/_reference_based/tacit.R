log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Annotate Xenium from Chromium using RCTD

library(Seurat)
library(Matrix)
library(dplyr)
library(spacexr)
library(arrow)
library(data.table)
library(future)

if (!requireNamespace("TACIT", quietly = TRUE)) {
  devtools::install_github("huynhkl953/TACIT@v1.0.0")
}
library(TACIT)

options(future.globals.maxSize = snakemake@params[["future_globals_maxSize"]])

# Load common reference-based parameters
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_header.R")

# parameters specific for RCTD and panel (or disease)
CELL_MIN_INSTANCE <- snakemake@params[["cell_min_instance"]]
ref_mode <- snakemake@params[["mode"]]
N_top <- snakemake@params[["N_binary_markers_per_cell_type"]]
N_PCs <- snakemake@params[["N_PCs"]]
r <- snakemake@params[["r"]]

###### end of snakemake params  ###### 


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'data' 
test_layer <- 'counts' 

# Load reference and query data
snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_load_data.R") # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods
rctd_path    <- snakemake@input[["rctd"]] 

expr <- GetAssayData(xe, assay = xe_assay, layer = test_layer) %>% t()
xe_chrom_common_genes <- intersect(rownames(xe), rownames(chrom))


if(ref_mode == "binary"){
  # Generate reference object
  snakemake@source("../../../scripts/_cell_type_annotation/_reference_based/_generate_reference_obj.R")
  ref <- generate_marker_binary_matrix(chrom[xe_chrom_common_genes,], celltypes = chrom[[annotation_level]] %>% pull, top_n = N_top, min_padj = 0.05) 
  message("Binary marker matrix has been succesfully generated. \n")
} else if(ref_mode == "continuous"){
  rctd <- readRDS(rctd_path)
  ref <- rctd@cell_type_info[[1]][[1]] %>% t() %>% as.matrix()
} else {
  stop(paste("Unknown `ref_mode`", ref_mode, ", available modes are `binary` and `continuous`"))
}
ref_df <- data.frame(cell_type = rownames(ref), ref)
expr <- expr[, colnames(ref)]

message("Starting TACIT... \n")
tct <- TACIT(data_expression = expr, Signature = ref_df, r=r, p=N_PCs)
message("TACIT Done ... \n")

ticit_weights <- tct[[2]][, -c(1:2)]
labels_df <- tct[[3]]
ticit_labels <- data.frame(cell_id = rownames(expr), label = labels_df[, "mem"], N_types = rowSums(ticit_weights))
#rownames(ticit_labels) <- ticit_labels$cell_id
#xe <- AddMetaData(xe, ticit_labels)

ticit_weights$cell_id <- rownames(expr)
ticit_weights <- ticit_weights %>% select(cell_id, everything())


# Save annotation
message("Writing results ... \n")
saveRDS(tct, snakemake@output[[1]])
write_parquet(ticit_labels, snakemake@output[[2]])
write_parquet(ticit_weights, snakemake@output[[3]])
message("Done Writing results ... \n")
