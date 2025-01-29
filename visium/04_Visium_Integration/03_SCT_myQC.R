library(Seurat)
library(ggplot2)
library(ggridges)
library(gtable)
library(patchwork)
library(dplyr)
library(stringr)

# disease = "breast"
disease = "lung"
# disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/16_Visium_Integration/Owkin_Script/00_params.R")


## Before SpotClean ------------------------------------------------------
datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
savepath.filtered = paste0(datapath.filtered, "/Results")
indication = paste0(save_names, "_qcd.rds")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-SCTnoSpotClean_filtered.pdf")
save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTnoSpotClean_filtered.rds")
datapath = datapath.filtered
savepath = savepath.filtered


# ## After SpotClean --------------------------------------------------------
# datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
# savepath.spotclean = paste0(datapath.spotclean, "/Results")
# indication = paste0("SpotClean_", save_names, "_qcd.rds")
# UMAP_name = paste0("/", str_to_title(disease), "-UMAP-SCTpostSpotClean.pdf")
# save_rds_name = paste0("/", str_to_title(disease), "-merge-SCTpostSpotClean.rds")
# datapath = datapath.spotclean
# savepath = savepath.spotclean


#########################################################################
indication.list=lapply(indication, function(filename){ 
  
  rds = filename
  samplename= gsub("SpotClean_","", gsub("_qcd.rds","", filename)) # works for both with or without "SpotClean_", would just ignore
  print(rds)
  if(file.exists(file.path(datapath, rds))){
    se_obj = readRDS(file.path(datapath, rds))
    # ## add hard QC
    # se_obj = subset(se_obj, nCount_Spatial >200) # Previously QCed with the same rule as SCE objects for decon and baye
    se_obj[["sample_id"]] = samplename
    return(se_obj)  
  }
  
})
names(indication.list) = save_names


#########################################################################
# SCtransform 
indication.list <- lapply(indication.list, FUN = function(x) {
  # x <- SCTransform(x, assay = "Spatial", verbose = FALSE, vst.flavor = "v1")#, do.scale = FALSE, do.center = FALSE)
  x <- SCTransform(x, assay = "Spatial", verbose = FALSE)#, do.scale = FALSE, do.center = FALSE) # no glmGamPoi, but still v2
})
indication.sct.features = SelectIntegrationFeatures(object.list = indication.list, nfeatures = 3000)
indication.sctmerge = merge(indication.list[[1]], y = indication.list[2:length(indication.list)], project = disease, merge.data=TRUE)
VariableFeatures(indication.sctmerge) <- indication.sct.features
indication.sctmerge <- RunPCA(object = indication.sctmerge, assay = "SCT", npcs = 50)
indication.sctmerge <- FindNeighbors(object = indication.sctmerge, assay = "SCT", reduction = "pca", dims = 1:50)
# indication.sctmerge <- FindClusters(object = indication.sctmerge, resolution = 0.4)
indication.sctmerge <- FindClusters(object = indication.sctmerge, resolution = 0.8)
indication.sctmerge=RunUMAP(object = indication.sctmerge, assay = "SCT", reduction = "pca", dims = 1:50)

pdf(paste0(savepath, UMAP_name), width = 7, height = 6)
DimPlot(indication.sctmerge, label=TRUE, group.by="Level4_decon_max") + labs(y= "UMAP_2", x = "UMAP_1")
DimPlot(indication.sctmerge, label=TRUE, group.by="Region") + labs(y= "UMAP_2", x = "UMAP_1")
print(DimPlot(indication.sctmerge, label=TRUE) + labs(y= "UMAP_2", x = "UMAP_1") )
print(DimPlot(indication.sctmerge, label=FALSE, group.by="sample_id") + labs(y= "UMAP_2", x = "UMAP_1"))
dev.off()

saveRDS(indication.sctmerge, file=paste0(savepath, save_rds_name))

# indication.sctmerge <- readRDS("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/breast/filtered/Results/Breast-merge-SCTnoSpotClean_filtered.rds")
