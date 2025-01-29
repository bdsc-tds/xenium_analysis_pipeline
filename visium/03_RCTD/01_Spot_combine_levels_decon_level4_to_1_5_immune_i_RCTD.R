library(dplyr)
library(tidyverse)
# Merge RCTD level 4 to level 1.5 immune
# decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/RCTD"
# decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/RCTD_DEgenes"
decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/C2L"
# decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/C2L_DEgenes"
# decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/Level4/AnchorIntegration"

# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD/level1_5_immune"
# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD_DEgenes/level1_5_immune"
agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L/level1_5_immune"
# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L_DEgenes/level1_5_immune"
# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/AnchorIntegration/level1_5_immune"

save_names_br <- c("B1_2",      "B1_4",      "B2_2",      "B3_2",      "B4_2"     )
save_names_lu <- c("L1_2",      "L1_4",      "L2_2",      "L3_2",      "L4_2"     )
save_names_dl <- c("DLBCL_1", "DLBCL_2", "DLBCL_3", "DLBCL_4", "DLBCL_5", "DLBCL_6")

# Breast -------------------------------------------------------------------
merge_Save <- function(name, decon_path, agg_decon_path, disease = "breast"){
  # test <- read.csv(file.path(decon_path, paste0(name, "_spot_level4_RCTD.csv"))) # RCTD*
  test <- read.csv(file.path(decon_path, paste0(name, ".csv"))) # C2L*
  # test_ <- data.frame(t(read.csv(file.path(decon_path, paste0(name, ".csv"))))) %>% janitor::row_to_names(row_number = 1); # AnchorIntegration
  # test <- data.frame(lapply(test_, as.numeric)) # AnchorIntegration
  # rownames(test) <- str_replace(rownames(test_), "\\.", "-") # AnchorIntegration
  # colnames(test) <- str_replace_all(colnames(test), "\\.", "_") # AnchorIntegration
  # test$X <- rownames(test) # AnchorIntegration
  # if(disease %in% c("breast", "lung")){ # AnchorIntegration
  #   test$Fibroblast <- test$Fibroblast + test$Fibroblast_B3; test$Fibroblast_B3 <- NULL
  # }

  if(disease == "breast"){
    # ## commented out for C2L all genes or AnchorIntegration
    # if(name == "B3_2"){
    # test <- test %>%
    #   dplyr::rename(Fibroblast = Fibroblast_B3)
    # }
    
    # C2L all genes -----
    if(name == "B3_2"){
      test <- test %>%
        dplyr::rename(Fibroblast = Fibroblast_B3)
    }else if(name %in% c("B1_2", "B1_4", "B2_2")){
      test <- test %>%
        mutate(Fibroblast = Fibroblast + Fibroblast_B3) %>%
        select(-Fibroblast_B3)
    }
    # -------------------

    test_merged <- test %>%
      rowwise() %>%
      mutate(`B cells` = sum(c_across(c(B_cell, B_plasma_IGHA1, B_plasma_IGHG1, B_plasma_IGHG3, B_plasma_IGHM, B_plasma_IGKC, B_plasma_IGLC1))),   
             Stroma = sum(across(c(Endothelia_vascular, Fibroblast, Muscle_smooth, Pericyte))),
             Tumor = sum(c_across(starts_with("Tu_B"))),
             `Myeloid else` = sum(c_across(c(DC_1, DC_2, DC_activated, DC_pc, Granulocyte, Mast_cell))),
             NK = NK,
             Macrophage = sum(c_across(c(Macrophage, Monocyte))),
             `T cells` = sum(c_across(c(T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg, TNK_dividing)))) %>%
      select(X, `B cells`, Stroma, Tumor, `Myeloid else`, NK, Macrophage, `T cells`)
    
  }else if(disease == "lung"){
    test_merged <- test %>%
      rowwise() %>%
      mutate(`B cells` = sum(c_across(c(B_cell, B_plasma_IGHA1, B_plasma_IGHG1, B_plasma_IGHG3, B_plasma_IGHM, B_plasma_IGKC, B_plasma_IGLC1))),   
             Stroma = sum(across(c(Endothelia_vascular, Fibroblast, Muscle_smooth, Pericyte))),
             Epithelia = Epi_lung,
             Tumor = sum(c_across(starts_with("Tu_L"))),
             `Myeloid else` = sum(c_across(c(DC_1, DC_2, DC_activated, DC_pc, Granulocyte, Mast_cell))),
             NK = NK,
             Macrophage = sum(c_across(c(Macrophage, Monocyte))),
             `T cells` = sum(c_across(c(T_CD4, T_CD8_exhausted, T_CTL, T_CXCL13, T_reg, TNK_dividing)))) %>%
      select(X, `B cells`, Stroma, Epithelia, Tumor, `Myeloid else`, NK, Macrophage, `T cells`)
    
  }else{ # DLBCL
    test_merged <- test %>%
      rowwise() %>%
      mutate(`B cells` = B_plasma,   
             Stroma = sum(across(c(Endothelia, Fibro_Muscle, any_of(c("Pericyte", "pericytes"))))),
             Epithelia = sum(c_across(starts_with("Epi_"))), 
             Tumor = sum(c_across(starts_with("Tu_D"))),
             `Myeloid else` = sum(c_across(c(DC_1, DC_2, DC_pc))),
             NK = NK,
             Macrophage = Mono_Macro,
             `T cells` = sum(c_across(c(T_CD4, T_CD4_reg, T_CD8, T_dividing)))) %>%
      select(X, `B cells`, Stroma, Epithelia, Tumor, `Myeloid else`, NK, Macrophage, `T cells`)
  }
  
  # write.csv(test_merged, file.path(agg_decon_path, paste0(name, "_spot_merged_RCTD.csv")), row.names = FALSE) # RCTD*
  write.csv(test_merged, file.path(agg_decon_path, paste0(name, ".csv")), row.names = FALSE) # C2L*
}

lapply(save_names_br, merge_Save, decon_path = decon_path, agg_decon_path = agg_decon_path, disease = "breast")
lapply(save_names_lu, merge_Save, decon_path = decon_path, agg_decon_path = agg_decon_path, disease = "lung")
lapply(save_names_dl, merge_Save, decon_path = decon_path, agg_decon_path = agg_decon_path, disease = "dlbcl")




