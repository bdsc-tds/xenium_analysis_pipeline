library(stringr)
library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)


# Helper ------------------------------------------------------------------
read_and_plot <- function(savepath, save_rds_name){
  seu <- readRDS(paste0(savepath, save_rds_name))
  p1 <- DimPlot(seu, label=TRUE, group.by="Region", pt.size = 0.85, label.size = 8) + 
    labs(y= "UMAP_2", x = "UMAP_1") + 
    theme(legend.position = "none", 
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) + 
    ggtitle("")
  p2 <- DimPlot(seu, label=TRUE, group.by="Level4_decon_max", pt.size = 0.85, label.size = 6) + 
    labs(y= "UMAP_2", x = "UMAP_1") + 
    theme(legend.position = "none", 
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.title = element_blank()) + 
    ggtitle("")
  p <- p1 + theme(plot.margin = margin(0, 50, 50, 0)) | p2
  return(p)
}
# p <- read_and_plot(savepath, save_rds_name)

###########################################################################
disease_list = c("breast", "lung", "dlbcl")
# disease = "breast"
# disease = "lung"
# disease = "dlbcl"

for(disease in disease_list){
  ## SCT #################################################################
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
p_combined <- p_breast_4 / p_lung_4 / p_dlbcl_4

figpath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Results/Manuscript_Figures_Final/SuppFig/Visium_Integration_UMAP_Heatmap"
pdf(file.path(figpath, "Visium_UMAP_Patho_Decon.pdf"), width = 26, height = 40)
print(p_combined)
dev.off()

