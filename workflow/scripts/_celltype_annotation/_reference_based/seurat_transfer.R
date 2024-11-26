# Annotate Xenium from Chromium using Seurat Transfer

library(Seurat)
library(dplyr)

# Load common reference-based parameters
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_header.R")) 


# Annotation method specific, *should not be modified*, so not in snakemake
ref_layer  <- 'data' 
test_layer <- 'data' 

# Load reference and query data
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_load_data.R")) # this does not look nice, but like this we do not duplicate code and making sure data are loaded and processed the same way for all the methods

## Make sure chrom data are log-normalized
DefaultAssay(chrom) <- ref.assay
chrom               <- NormalizeData(chrom)

# Generate reference object
source(snakemake@source("../../../scripts/_celltype_annotation/_reference_based/_generate_reference_obj.R")) 

# ? Should chrom data be re-normalized besed on shared genes???

ref_labels <- chrom@meta.data %>% pull(annotation_level) 

xe_chrom_common_genes <- intersect(rownames(xe), rownames(chrom))

# Annotate with Seurat transfer 

dims <- snakemake@params[["dims"]] %||% 1:50

anchors       <- FindTransferAnchors(reference = chrom, query = xe, features = xe_chrom_common_genes, dims = dims)
predictions   <- TransferData(anchorset = anchors, refdata = ref_labels, dims = dims)
xe            <- AddMetaData(xe, metadata = predictions)


scores                <- predictions %>% select(!c("predicted.id", "prediction.score.max"))
colnames(scores)      <- gsub("prediction.score.", "", colnames(scores))
labels                <- predictions$predicted.id
rownames(scores)      <- colnames(xe)
names(labels)         <- colnames(xe)

# Save annotation
saveRDS(predictions, snakemake@output[[1]])
write.csv(labels, snakemake@output[[2]])
write.csv(scores, snakemake@output[[3]]) 

