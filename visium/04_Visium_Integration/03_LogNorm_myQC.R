library(Seurat)
library(ggplot2)
library(ggridges)
library(gtable)
library(patchwork)
library(dplyr)
library(stringr)

# disease = "breast"
# disease = "lung"
disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/16_Visium_Integration/Owkin_Script/00_params.R")


# ## Before SpotClean ------------------------------------------------------
# datapath.filtered = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/filtered")
# savepath.filtered = paste0(datapath.filtered, "/Results")
# indication = paste0(save_names, "_qcd.rds")
# UMAP_name = paste0("/", str_to_title(disease), "-UMAP-logNormnoSpotClean_filtered.pdf")
# save_rds_name = paste0("/", str_to_title(disease), "-merge-logNormnoSpotClean_filtered.rds")
# datapath = datapath.filtered
# savepath = savepath.filtered


## After SpotClean --------------------------------------------------------
datapath.spotclean = paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean")
savepath.spotclean = paste0(datapath.spotclean, "/Results")
indication = paste0("SpotClean_", save_names, "_qcd.rds")
UMAP_name = paste0("/", str_to_title(disease), "-UMAP-logNormpostSpotClean.pdf")
save_rds_name = paste0("/", str_to_title(disease), "-merge-logNormpostSpotClean.rds")
datapath = datapath.spotclean
savepath = savepath.spotclean


#########################################################################
indication.list=lapply(indication, function(filename){ 
  
  rds = filename
  samplename= gsub("SpotClean_","", gsub("_qcd.rds","", filename)) # works for both with or without "SpotClean_", would just ignore
  print(rds)
  if(file.exists(file.path(datapath, rds))){
    se_obj = readRDS(file.path(datapath, rds))
    # ## add hard QC
    # se_obj = subset(se_obj, nCount_Spatial >200)
    se_obj[["sample_id"]] = samplename
    return(se_obj)  
  }
  
})
names(indication.list) = save_names


#########################################################################
# LogNorm 
indication.list <- lapply(indication.list, FUN = function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)  
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})
indication.logn.features = SelectIntegrationFeatures(object.list = indication.list, nfeatures = 3000)
indication.lognmerge = merge(indication.list[[1]], y = indication.list[2:length(indication.list)], project = disease, merge.data=TRUE)
VariableFeatures(indication.lognmerge) <- indication.logn.features
indication.lognmerge <- ScaleData(indication.lognmerge)
indication.lognmerge <- RunPCA(object = indication.lognmerge, assay = "Spatial", npcs = 50)
indication.lognmerge <- FindNeighbors(object = indication.lognmerge, assay = "Spatial", reduction = "pca", dims = 1:50)
indication.lognmerge <- FindClusters(object = indication.lognmerge, resolution = 0.8)
indication.lognmerge=RunUMAP(object = indication.lognmerge, assay = "Spatial", reduction = "pca", dims = 1:50)

pdf(paste0(savepath, UMAP_name), width = 7, height = 6)
DimPlot(indication.lognmerge, label=TRUE, group.by="Level4_decon_max") + labs(y= "UMAP_2", x = "UMAP_1")
DimPlot(indication.lognmerge, label=TRUE, group.by="Region") + labs(y= "UMAP_2", x = "UMAP_1")
print(DimPlot(indication.lognmerge, label=TRUE))
print(DimPlot(indication.lognmerge, label=TRUE, group.by="sample_id"))
dev.off()

saveRDS(indication.lognmerge, file=paste0(savepath, save_rds_name))


