log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(Seurat)
library(dplyr)
library(stringr)
library(reader)
library(arrow)

# Parse command-line arguments - not sure if this section is correct, leaving for Senbai
args <- commandArgs(trailingOnly = TRUE)
seurat_files <- strsplit(args[1], " ")[[1]]  # Multiple Seurat files
annotation_file <- args[2]
output_file <- args[3]

# Load annotation data
annotations <- read_parquet(annotation_file)

# Function to process each Seurat object
process_seurat_metadata <- function(seurat_path) {
    path_parts <- unlist(strsplit(seurat_path, "/"))
    # this should be taken from the annotation file instead!
    segmentation <- path_parts[1]
    condition <- path_parts[2]
    panel <- path_parts[3]  # Extract panel
    donor <- path_parts[4]
    sample_id <- path_parts[5]
    reference <- path_parts[8]
    annotation_level <- path_parts[10]

    seurat_obj <- readRDS(seurat_path)

    # Add metadata
    seurat_obj <- AddMetaData(seurat_obj, metadata = annotations)

    # Extract metadata
    metadata <- seurat_obj@meta.data

    # Summarize count per sample
    sample_summary <- metadata %>%
        summarise(
            sample_id = unique(sample_id),
            panel = panel,
            segmentation = segmentation,
            reference = reference,
            doublet_detection_method = "scDblFinder",
            annotation_level = annotation_level,
            Count = n()
        )
    
    return(sample_summary)
}

# Process all Seurat objects and combine results
all_summaries <- lapply(seurat_files, process_seurat_metadata)
final_summary <- bind_rows(all_summaries)

# Group by panel and write separate output files
final_summary %>%
    group_by(panel) %>%
    write.csv(output_file, row.names = FALSE)
