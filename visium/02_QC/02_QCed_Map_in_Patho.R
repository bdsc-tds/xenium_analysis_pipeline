disease_list = c("lung", "breast", "dlbcl")

for(d in 1:3){
  disease = disease_list[d]
  source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
  
  # PCA ---------------------------------------------------------------------
  for(i in 1:nsamples){
    print(save_names[i])
    sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
    
    # merge patho
    if(save_names[i] != "B1_2"){
      vis_anno <- read.csv(file.path(patho_annopath, paste0(save_names[i], ".csv")))
      CD <- data.frame(colData(sce))
      CD$Barcode <- rownames(CD)
      
      CD <- CD %>%
        left_join(vis_anno, by = "Barcode")
      rownames(CD) <- CD$Barcode
      
      all(colnames(sce) == rownames(CD))
      
      colData(sce) <- as(CD, "DFrame")
    }else{ #  if(save_names[i] == c("B1_2"))
      sce$Region <- FALSE
      sce$Barcode <- colnames(sce)
    }
    
    saveRDS(sce, file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
  }
}