library(Seurat)
library(arrow)
source("workflow/scripts/utils/readwrite.R")

args <- commandArgs(trailingOnly = TRUE)
input_rds <- args[1]
output_dir <- args[2]

# Load Seurat object
seurat_obj <- readRDS(input_rds)

# Create output directory
dir.create(output_dir, showWarnings = FALSE)

### **1️⃣ Save Metadata**
metadata_path <- file.path(output_dir, "metadata.parquet")
metadata <- as.data.frame(seurat_obj@meta.data)
metadata$cell_id <- rownames(metadata) 
write_parquet(metadata, metadata_path)
message("✅ Saved metadata to: ", metadata_path)

### **2️⃣ Save Assays (counts, data, scale.data)**
for (assay_name in Assays(seurat_obj)) {
    message("🔹 Processing assay: ", assay_name)
    
    DefaultAssay(seurat_obj) <- assay_name  # Set active assay

    for (slot in c("counts", "data", "scale.data")) {

        message("  - Extracting ", slot, " from ", assay_name)

        mat <- GetAssayData(seurat_obj, layer = slot, assay = assay_name)
        
        if((nrow(mat) > 0 && ncol(mat) > 0)){
            slot_parquet_path <- file.path(output_dir, paste0(assay_name, "_", slot, ".h5"))
            write10xCounts(path=slot_parquet_path, x=mat,type='HDF5')
            message("  ✅ Saved: ", slot_parquet_path)

        } else{
            message("  ❌ Skipping ", slot, " (empty matrix in ", assay_name, ")")

        }

    }
}

### **3️⃣ Save Dimensional Reductions (UMAP, PCA, etc.)**
for (red_name in Reductions(seurat_obj)) {
    message("🔹 Processing reduction: ", red_name)
    red_obj <- seurat_obj[[red_name]]  # Extract reduction
    
    # Save embeddings
    if ("cell.embeddings" %in% slotNames(red_obj)) {
        emb_path <- file.path(output_dir, paste0(red_name, "_embeddings.parquet"))
        write_parquet(as.data.frame(red_obj@cell.embeddings), emb_path)
        message("  ✅ Saved embeddings: ", emb_path)
    }
    
    # Save feature loadings (if available)
    # if ("feature.loadings" %in% slotNames(red_obj)) {
    #     load_path <- file.path(output_dir, paste0(red_name, "_loadings.parquet"))
    #     write_parquet(as.data.frame(red_obj@feature.loadings), load_path)
    #     message("  ✅ Saved feature loadings: ", load_path)
    # }
    
    # Save standard deviation (if available)
    # if ("stdev" %in% slotNames(red_obj)) {
    #     stdev_path <- file.path(output_dir, paste0(red_name, "_stdev.parquet"))
    #     write_parquet(as.data.frame(red_obj@stdev), stdev_path)
    #     message("  ✅ Saved stdev: ", stdev_path)
    # }
}

file.create(file.path(output_dir,".done"))
message("🎉 All data have been exported!")
