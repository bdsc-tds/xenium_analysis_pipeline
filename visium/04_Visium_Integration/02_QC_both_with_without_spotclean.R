library(SingleCellExperiment)
library(SpatialExperiment)
library(dplyr)
# Map over QCed rownames and colnames from filtered (counts) SCE QC criteria in 
# "/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/02_QC_Check.R", same criteria as for Decon and Baye

# Breast ----------------------------------------------------
disease <- "breast"
sample_name <- c("B4_2", "B3_2", "B2_2", "B1_2", "B1_4")
# dims <- list(c(18085, 1850), c(18085, 2220), c(18085, 2569), c(18085, 1897), c(18085, 1999))
nsamples = 5

# Lung ------------------------------------------------------
disease <- "lung"
sample_name <- c("L4_2", "L3_2", "L2_2", "L1_2", "L1_4")
# dims <- list(c(18085, 2958), c(18085, 2443), c(18085, 842), c(18085, 874), c(18085, 944))
nsamples = 5

# DLBCL ----------------------------------------------------
disease <- "dlbcl"
sample_name <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
# dims <- list(c(18085, 4710), c(18085, 4121), c(18085, 4951), c(18085, 1457), c(18085, 1782), c(18085, 1741))
nsamples = 6


qcd_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd/")
# Read both filtered and spotcleaned object
integration_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/"
filtered_path <- paste0(integration_path, disease, "/filtered/")
spcl_path <- paste0(integration_path, disease, "/spotclean/")

for(n in 1:nsamples){
  samplename <-sample_name[n]
  
  qcd_sce <- readRDS(file.path(qcd_path, paste0(samplename, "_qcd.rds")))
  qcd_sce_df <- data.frame(colData(qcd_sce)) %>%
    select(Barcode, Region, Level4_decon_max)
  
  
  # Map in pathology and decon -------------------------------------------
  filtered_seu <- readRDS(file.path(filtered_path, paste0(samplename, ".rds")))
  filtered_seu_qcd <- filtered_seu[rownames(filtered_seu) %in% rownames(qcd_sce), colnames(filtered_seu) %in% colnames(qcd_sce)]
  filtered_seu_qcd_df <- data.frame(Barcode = colnames(filtered_seu_qcd)) %>%
    left_join(qcd_sce_df)
  
  filtered_seu_qcd$Region <- filtered_seu_qcd_df$Region
  filtered_seu_qcd$Level4_decon_max <- filtered_seu_qcd_df$Level4_decon_max

  saveRDS(filtered_seu_qcd, file.path(filtered_path, paste0(samplename, "_qcd.rds")))
  
  
  # Map in pathology and decon -------------------------------------------
  spcl_seu <- readRDS(file.path(spcl_path, paste0("SpotClean_", samplename, ".rds")))
  spcl_seu_qcd <- spcl_seu[rownames(spcl_seu) %in% rownames(qcd_sce), colnames(spcl_seu) %in% colnames(qcd_sce)]
  spcl_seu_qcd_df <- data.frame(Barcode = colnames(spcl_seu_qcd)) %>%
    left_join(qcd_sce_df)
  
  spcl_seu_qcd$Region <- spcl_seu_qcd_df$Region
  spcl_seu_qcd$Level4_decon_max <- spcl_seu_qcd_df$Level4_decon_max
  
  saveRDS(spcl_seu_qcd, file.path(spcl_path, paste0("SpotClean_", samplename, "_qcd.rds")))
}




