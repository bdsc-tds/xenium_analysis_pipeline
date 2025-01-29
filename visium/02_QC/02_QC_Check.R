# disease = "breast"
disease = "lung"
# disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")

# PCA ---------------------------------------------------------------------
i = 4
for(i in 1:nsamples){
  sce <- readRDS(file.path(save_path, paste0(sample_names[i], ".rds")))
  dim(sce)
  
  ### Raw - QC should be on raw
  print(save_names[i])
  
  ## Spot-wise
  sce <- addPerCellQC(sce, subsets=list(Mito=grepl("^MT-", rownames(sce))))
  sce$subsets_Mito_percent <- ifelse(is.na(sce$subsets_Mito_percent), 0, sce$subsets_Mito_percent)
  
  # Mito
  tail(sort(sce$subsets_Mito_percent), 20)
  sce$mito_drop <- sce$subsets_Mito_percent > 22
  plotSpotQC(sce, type = "spots", annotate = "mito_drop", 
             x_coord = "array_col", y_coord = "array_row")
  
  # Sum
  head(sort(sce$sum), 50)
  sce$libsize_drop <- sce$sum < 100
  plotSpotQC(sce, type = "spots", annotate = "libsize_drop", 
             x_coord = "array_col", y_coord = "array_row")
  
  # Patho
  if(save_names[i] != "B1_2"){
    vis_anno <- read.csv(file.path("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium/Annotations", paste0(save_names[i], ".csv")))
    sce$patho_exclude <- colnames(sce) %in% vis_anno$Barcode[vis_anno$Region == "Exclude"] # Merge in patho Region and | is.na(sce$Region)
  }else{ #  if(save_names[i] == c("B1_2"))
    sce$patho_exclude <- FALSE
  }

  plotSpotQC(sce, type = "spots", annotate = "patho_exclude", 
             x_coord = "array_col", y_coord = "array_row")
  
  # Edge location
  range(sce$array_col)
  range(sce$array_row)
  sort(sce$array_row[sce$array_col == 116])
  sort(sce$array_row[sce$array_col == 24])
  sce$array_col[sce$array_row == 24]
  
  # sce$edge <- (sce$array_col == 0 & sce$array_row == 22) |   # B1_2
  #   (sce$array_col == 1 & sce$array_row == 21)
  # sce$edge <- (sce$array_col == 125 & sce$array_row == 59) | # B1_4
  #   (sce$array_col == 125 & sce$array_row == 71)
  # sce$edge <- FALSE                                          # B2_2
  # sce$edge <- (sce$array_col == 126 & sce$array_row == 68)   # B3_2
  # sce$edge <- FALSE                                          # B4_2
  # sce$edge <- FALSE                                          # L1_2
  # sce$edge <- FALSE                                          # L1_4
  # sce$edge <- (sce$array_col == 48 & sce$array_row == 74)    # L2_2
  # sce$edge <- FALSE                                          # L3_2
  # sce$edge <- FALSE                                          # L4_2
  # sce$edge <- FALSE                                          # DLBCL_1
  # sce$edge <- FALSE                                          # DLBCL_2
  # sce$edge <- FALSE                                          # DLBCL_3
  # sce$edge <- FALSE                                          # DLBCL_4
  # sce$edge <- (sce$array_col == 116 & sce$array_row == 72) | # DLBCL_5
  #   (sce$array_col == 120 & sce$array_row == 64) | 
  #   (sce$array_col == 123 & sce$array_row == 75)
  # sce$edge <- FALSE                                          # DLBCL_6
  plotSpotQC(sce, type = "spots", annotate = "edge", 
             x_coord = "array_col", y_coord = "array_row")

  ## Gene-wise
  sce <- addPerFeatureQC(sce)
  rowData(sce)$low_abungene <- rowSums(assays(sce)$counts > 0) < 20
  print(paste0(length(which(rowData(sce)$low_abungene)), " genes")) 
  
  ## Plot
  sce$all_drop <- (sce$mito_drop | sce$libsize_drop | sce$edge | sce$patho_exclude) 
  plotSpotQC(sce, type = "spots", annotate = "all_drop", 
             x_coord = "array_col", y_coord = "array_row")
  
  plotFeatureQC(sce, type = "violin", metric_x = "detected", annotate = "low_abungene")
  
  ## Drop
  sce <- sce[!rowData(sce)$low_abungene, !sce$all_drop]
  dim(sce)
  
  save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_qcd/", disease, "_qcd")
  saveRDS(sce, file.path(save_bs_path, paste0(save_names[i], "_qcd.rds")))
  
}

## sample dim 
# L1_2: 13682   870
# L1_4: 14866   944
# L2_2: 11795   822
# L3_2: 17919  2378
# L4_2: 17863  2754

# B1_2: 17293  1893
# B1_4: 16668  1988
# B2_2: 17626  2560
# B3_2: 17883  2156
# B4_2: 13350  1652

# DLBCL_1: 18020  3739
# DLBCL_2: 17993  4114
# DLBCL_3: 18000  2019
# DLBCL_4: 18007  1436
# DLBCL_5: 18011  1733
# DLBCL_6: 18008  1650
