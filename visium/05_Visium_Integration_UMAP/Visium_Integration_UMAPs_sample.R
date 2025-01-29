library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)


# Helper ------------------------------------------------------------------
read_and_plot <- function(savepath, save_rds_name){
  seu <- readRDS(paste0(savepath, save_rds_name))
  p <- DimPlot(seu, label=TRUE, group.by="sample_id", label.size = 8) + 
    labs(y= "UMAP_2", x = "UMAP_1") + ggtitle("") + 
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank(),
          legend.text = element_text(size=15)) 
  return(p)
}
# p <- read_and_plot(savepath, save_rds_name)

###########################################################################
disease_list = c("breast", "lung", "dlbcl")
# disease = "breast"
# disease = "lung"
# disease = "dlbcl"

for(disease in disease_list){
  ## LogNorm ##############################################################
  ## Before SpotClean -----------------------------------------------------
  datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
  savepath.filtered = paste0(datapath.filtered, "/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-logNormnoSpotClean_filtered.rds")
  savepath = savepath.filtered
  pipeline = 1
  
  p <- read_and_plot(savepath, save_rds_name)
  assign(paste0("p_", disease, "_", pipeline), p)
  
  ## After SpotClean ------------------------------------------------------
  datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
  savepath.spotclean = paste0(datapath.spotclean, "/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-logNormpostSpotClean.rds")
  savepath = savepath.spotclean
  pipeline = 2
  
  p <- read_and_plot(savepath, save_rds_name)
  assign(paste0("p_", disease, "_", pipeline), p)
  
  ## SCT #################################################################
  # Before spotclean -----------------------------------------------------
  datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
  savepath.filtered = paste0(datapath.filtered, "/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTnoSpotClean_filtered.rds")
  savepath = savepath.filtered
  pipeline = 3
  
  p <- read_and_plot(savepath, save_rds_name)
  assign(paste0("p_", disease, "_", pipeline), p)
  
  # After spotclean ------------------------------------------------------
  datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
  savepath.spotclean = paste0(datapath.spotclean, "/Results")
  save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
  savepath = savepath.spotclean
  pipeline = 4
  
  p <- read_and_plot(savepath, save_rds_name)
  assign(paste0("p_", disease, "_", pipeline), p)
}

#########################################################################
p_combined <- (p_breast_1 + theme(plot.margin = margin(0, 50, 50, 0)) | p_lung_1 + theme(plot.margin = margin(0, 50, 50, 50)) | p_dlbcl_1 + theme(plot.margin = margin(0, 0, 50, 50))) /
  (p_breast_2 + theme(plot.margin = margin(0, 50, 50, 0)) | p_lung_2 + theme(plot.margin = margin(0, 50, 50, 50)) | p_dlbcl_2 + theme(plot.margin = margin(0, 0, 50, 50))) /
  (p_breast_3 + theme(plot.margin = margin(0, 50, 50, 0)) | p_lung_3 + theme(plot.margin = margin(0, 50, 50, 50)) | p_dlbcl_3 + theme(plot.margin = margin(0, 0, 50, 50))) /
  (p_breast_4 + theme(plot.margin = margin(0, 50, 0, 0)) | p_lung_4 + theme(plot.margin = margin(0, 50, 0, 50)) | p_dlbcl_4 + theme(plot.margin = margin(0, 0, 0, 50))) 

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_UMAP_Heatmap"
pdf(file.path(figpath, "Visium_UMAP_Samples.pdf"), width = 26, height = 30)
print(p_combined)
dev.off()
