library(Seurat)
library(SpotClean)
library(SummarizedExperiment)
library(cowplot)

# Breast ----------------------------------------------------
disease <- "breast"
sample_id <- c("1FHZ", "1GVR", "1256", "OPHI", "OPHI")
section_id <- c("V10_B4_2_1FHZ", "V8_B3_2_1GVR", "V5_B2_2_1256", "V4_B1_2_OPHI", "V6_B1_4_OPHI")
sample_name <- c("B4_2", "B3_2", "B2_2", "B1_2", "B1_4")
# dims <- list(c(18085, 1850), c(18085, 2220), c(18085, 2569), c(18085, 1897), c(18085, 1999))
nsamples = 5

# Lung ------------------------------------------------------
disease <- "lung"
sample_id <- c("1GA2", "1G73", "0WMU", "0PSV", "0PSV")
section_id <- c("V13_L4_2_1GA2", "V12_L3_2_1G73", "V11_L2_2_0WMU", "V7_L1_2_0PSV", "V9_L1_4_0PSV")
sample_name <- c("L4_2", "L3_2", "L2_2", "L1_2", "L1_4")
# dims <- list(c(18085, 2958), c(18085, 2443), c(18085, 842), c(18085, 874), c(18085, 944))
nsamples = 5

# DLBCL ----------------------------------------------------
disease <- "dlbcl"
sample_id <- c("V14_DLBCL_1", "V15_DLBCL_2", "V16_DLBCL_3", "V17_DLBCL_4", "V18_DLBCL_5", "V19_DLBCL_6")
sample_name <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")
# dims <- list(c(18085, 4710), c(18085, 4121), c(18085, 4951), c(18085, 1457), c(18085, 1782), c(18085, 1741))
nsamples = 6

datapath.filtered <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered/")

for(n in 1:nsamples){
  if(disease %in% c("breast", "lung")){
    data_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot/",
                        disease, "/",
                        sample_id[n],
                        "/visium/",
                        section_id[n],
                        "/outs/")
  }else if(disease == "dlbcl"){
    if(n %in% c(1,3)){
      data_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/dlbcl_raw/Visium_D1_D3_hair/",
                          paste0("D", n), "/",
                          "outs/")
    }else{
      data_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/owkin_pilot/mosaic_pilot_2/visium/", 
                          sample_id[n], "/",
                          sample_id[n], "/",
                          "outs/")
    }
  }
  samplename <-sample_name[n]
  
  # filtered ----------------------------------------------------------------
  if(disease == "dlbcl" && n %in% c(1,3)){
    # Load Visium as Seurat object
    visium_Seurat <- Seurat::Load10X_Spatial(data.dir = data_path,                              # to get the filtered gene list 18085
                                             filename = "filtered_feature_bc_matrix.h5",
                                             image = Read10X_Image(paste0(data_path, "spatial/")))
    spe <- SpatialExperiment::read10xVisium(samples = data_path, data = "raw")                  # to get the in tissue barcodes 4710 for D1; 4951 for D3 (from pathologist)
    
    visium_Seurat_raw <- Seurat::Load10X_Spatial(data.dir = data_path,
                                                 filename = "raw_feature_bc_matrix.h5",         # to get the raw count matrix
                                                 image = Read10X_Image(paste0(data_path, "spatial/")))
    
    seu <- visium_Seurat_raw[rownames(visium_Seurat_raw) %in% rownames(visium_Seurat), colnames(visium_Seurat_raw) %in% colnames(spe)[spe$in_tissue]] # 18085  4710 Visium seurat filtered
    
  }else{
    seu <- Seurat::Load10X_Spatial(data.dir = data_path,
                                   filename = "filtered_feature_bc_matrix.h5",
                                   image = Read10X_Image(paste0(data_path, "spatial/")))
  }
  
  # SpatialFeaturePlot(seu, features = "nCount_Spatial") + theme(legend.position = "right")
  
  saveRDS(seu, paste0(datapath.filtered, samplename, ".rds"))
  
}

## In the R sessionInfo() container given by Sabrina, missing hdf5r package
# Error in Read10X_h5(filename = file.path(data.dir, filename), ...) : 
#   Please install hdf5r to read HDF5 files

## Therefore, copied the result of filtered object from my previous run with Seurat 4 + SpotClean 1.4.1
for(n in 1:6){
  samplename <-sample_name[n]
  print(dim(readRDS(paste0(datapath.filtered, samplename, ".rds"))))
}

## DLBCL correct dimension
# [1] 18085  4710
# [1] 18085  4121
# [1] 18085  4951
# [1] 18085  1457
# [1] 18085  1782
# [1] 18085  1741


