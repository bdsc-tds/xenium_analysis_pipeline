log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Load required libraries
library(tidyverse)
library(Seurat)

###### snakemake params  ######

query_path <- snakemake@input[[1]] # combo meta.data obj (will contain all Seurat metadata, including spot_class)

###### end of snakemake params  ###### 

# Function to process a segmentation dataset 
process_segmentation <- function(metadata, 
                                 count_col = "nCount", 
                                 feature_col = "nFeature",
                                 gene_panel = split_sample_id[2], # is this a function I should load first?
                                 donor_col = split_sample_id[3],
                                 sample_col = split_sample_id[4],
                                 sample_id = "sample_id",
                                 segmentation_col = "segmentation_id") {
  
  # Extract unique tissue and segmentation levels
  tissues <- unique(metadata[[sample_col]])
  segmentations <- unique(metadata[[segmentation_col]])
  
  # Initialize an empty list to store results
  all_results <- list()
  
  # Iterate over each tissue and segmentation method
  for (tissue in tissues) {
    for (segmentation in segmentations) {
      
      # Subset metadata object for the current tissue and segmentation method
      md_subset <- subset(metadata, subset = tissue == tissue & segmentation == segmentation)
      
      # Ensure the subset is not empty
      if (nrow(md_subset) > 0) {
        
        # Extract the count and feature data
        count_data <- md_subset[[count_col]]
        feature_data <- md_subset[[feature_col]]
        
        # Convert to numeric vectors (remove row names) before computing the mean
        count_data_numeric <- as.numeric(count_data)
        feature_data_numeric <- as.numeric(feature_data)
        
        # Calculate QC metrics for the current tissue and segmentation method
        number_of_cells <- ncol(md_subset)
        percentage_singlets <- ncol(md_subset %>% filter(spot_class == "singlet")) / number_of_cells * 100
        percentage_rejects <- ncol(md_subset %>% filter(spot_class == "reject")) / number_of_cells * 100
        
        # Assuming count_data_numeric and feature_data_numeric are numeric vectors, calculate the mean
        mean_nCount <- mean(count_data_numeric, na.rm = TRUE)
        mean_nFeature <- mean(feature_data_numeric, na.rm = TRUE)

        # Assuming count_data_numeric and feature_data_numeric are numeric vectors, calculate the median
        median_nCount <- median(count_data_numeric, na.rm = TRUE)
        median_nFeature <- median(feature_data_numeric, na.rm = TRUE)
        
        # Store the results in a list
        result <- data.frame(
          sample = sample_id, 
          segmentation = segmentation,
          donor = split_sample_id[3],
          number_of_cells = number_of_cells,
          mean_nCount = mean_nCount, 
          median_nCount = median_nCount,
          mean_nFeature = mean_nFeature,
          median_nFeature = median_nFeature,
          percentage_singlets = percentage_singlets, 
          percentage_rejects = percentage_rejects 
        )
        
        # Append the result to the list
        all_results[[length(all_results) + 1]] <- result
      }
    }
  }
  
  # Combine all the results into a single data frame
  final_data <- do.call(rbind, all_results)
}


# Running the function
metadata <- read.csv(query_path)
all_data <- process_segmentation(metadata)
write.csv(all_data, snakemake@output[[1]], row.names = FALSE)
