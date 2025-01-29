library(dplyr)
library(stringr)

# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD/level1_5_immune"
# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD_DEgenes/level1_5_immune"
agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L/level1_5_immune"
# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L_DEgenes/level1_5_immune"
# agg_decon_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/AnchorIntegration/level1_5_immune"

# agg_decon_long_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD/level1_5_immune_long"
# agg_decon_long_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/RCTD_DEgenes/level1_5_immune_long"
agg_decon_long_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L/level1_5_immune_long"
# agg_decon_long_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/C2L_DEgenes/level1_5_immune_long"
# agg_decon_long_path <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Visium_Decon/AggLevel4/AnchorIntegration/level1_5_immune_long"

get_Merged_Decon <- function(){
  decon_patho_gather <- NULL
  for(i in 1:nsamples){
    sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))
    patho <- data.frame(Region = sce$Region,
                        Barcode = sce$Barcode)
    
    # results <- read.csv(file.path(agg_decon_path, paste0(save_names[i], "_spot_merged_RCTD.csv"))) # RCTD*
    results <- read.csv(file.path(agg_decon_path, paste0(save_names[i], ".csv"))) # C2L*
    results <- results %>% dplyr::rename(Barcode = X)
    
    decon_patho_i <- results %>%
      left_join(patho, by = "Barcode")
    
    decon_patho_gather_i <- gather(decon_patho_i, key = "CellType", value = "Fraction", -c("Barcode", "Region"))
    decon_patho_gather_i$Section <- save_names[i]
    
    
    decon_patho_gather <- rbind(decon_patho_gather, decon_patho_gather_i)
  }
  
  decon_patho_gather$CellType <- str_replace(decon_patho_gather$CellType, "\\.", " ")
  return(decon_patho_gather)
}


# -------------------------------------------------------------------------
disease = "breast"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
decon_patho_gather_breast <- get_Merged_Decon()
write.csv(decon_patho_gather_breast, file.path(agg_decon_long_path, "vis_breast_decon_immune_long.csv"), row.names = FALSE)


# -------------------------------------------------------------------------
disease = "lung"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
decon_patho_gather_lung <- get_Merged_Decon()
write.csv(decon_patho_gather_lung, file.path(agg_decon_long_path, "vis_lung_decon_immune_long.csv"), row.names = FALSE)


# -------------------------------------------------------------------------
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))
decon_patho_gather_dlbcl <- get_Merged_Decon()
decon_patho_gather_dlbcl$Section <- case_when(decon_patho_gather_dlbcl$Section == "DLBCL_1" ~ "D1",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_2" ~ "D2",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_3" ~ "D3",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_4" ~ "D4",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_5" ~ "D5",
                                              decon_patho_gather_dlbcl$Section == "DLBCL_6" ~ "D6")
write.csv(decon_patho_gather_dlbcl, file.path(agg_decon_long_path, "vis_dlbcl_decon_immune_long.csv"), row.names = FALSE)


