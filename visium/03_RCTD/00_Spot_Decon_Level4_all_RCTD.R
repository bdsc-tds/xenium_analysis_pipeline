library(SpatialExperiment)
library(spacexr)

disease = "breast"
# disease = "lung"
# disease = "dlbcl"
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/Visium/Manuscript/01_params.R")
source("/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/ydong/Owkin_Pilot/Code/color_palette.R")

foldername <- ifelse(disease == "dlbcl", "DLBCL", str_to_title(disease))

## Get chromium by indication ------------------------------------------
chrompath <- "/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Chromium/Breast_Lung/For_manuscript_decon/" # breast / # lung
save_bs_path <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Intermediate/Visium_BayesSpace_raw/")

vis_rawmat <- "counts"

i = 1
# for(i in 1:nsamples){
# Chrom QCed --------------------------------------------------------------
if(disease == "breast"){
  chrom <- readRDS(file.path(chrompath, "chrom_breast_add_lung_healthy.rds"))
  tumor_classes <- c("Tu_B1_MUCL1", "Tu_B1_MUCL1_necrosis", "Tu_B1_MUCL1_transcription", "Tu_B3_CYP4F8", "Tu_B4_RHOB")
}else if(disease == "lung"){
  chrom <- readRDS(file.path(chrompath, "chrom_lung_add_breast_healthy.rds"))
  tumor_classes <- c("Tu_L1_SFTPB",                  "Tu_L2_FXYD2",        "Tu_L3_G0S2",           "Tu_L3_G0S2_immune_signature", 
                     "Tu_L4_KRT17_immune_signature", "Tu_L4_KRT17_mucous", "Tu_L4_KRT17_necrosis", "Tu_L4_KRT17_neutrophil_signature")
}else{ # DLBCL
  chrom <- readRDS(file.path(chrompath, "chrom_dlbcl.rds"))
  tumor_classes <- c("Tu_D1_LMO2",   "Tu_D1_RGS13",    "Tu_D1_SMIM14", 
                     "Tu_D2_mito",
                     "Tu_D3_BCL2A1", "Tu_D3_dividing", "Tu_D3_FAM3C", "Tu_D3_IGHD",
                     "Tu_D4_BCL7A",  "Tu_D4_PNN",      "Tu_D5_CCL22", "Tu_D6_BCL2")
}

if(disease %in% c("breast", "lung")){ # (btw B2_2 does not have tumor class)
  chrom_other_pt_tumor <- tumor_classes[!grepl(substr(save_names[i], 1, 2), 
                                               tumor_classes)]
  if(save_names[i] == "B3_2"){
    fibro_remove = "Fibroblast"
  }else{ # all other samples
    fibro_remove = "Fibroblast_B3"
  }
  
  if(save_names[i] == "B2_2"){ # Pull other breast tumor class into tumor
    chrom$Harmonised_Level4 <- ifelse(chrom$Harmonised_Level4 %in% tumor_classes, "Tu_B2", chrom$Harmonised_Level4)
  }
  to_remove = c(chrom_other_pt_tumor, fibro_remove)
}else{ # "dlbcl"
  chrom_other_pt_tumor <- tumor_classes[!grepl(paste0(substr(save_names[i], 1, 1), substr(save_names[i], 7, 7)), 
                                               tumor_classes)]
  to_remove = chrom_other_pt_tumor
}

chrom <- chrom[, !(chrom$Harmonised_Level4 %in% to_remove)]
chrom$Harmonised_Level4 <- droplevels(as.factor(chrom$Harmonised_Level4))
table(chrom$Harmonised_Level4)

cell_types <- as.factor(chrom$Harmonised_Level4); names(cell_types) <- colnames(chrom)
chrom_mat <- chrom@assays$RNA@counts 

## Get Spot gene expr by sample -----------------------------------------
sce <- readRDS(file.path(save_path_qcd, paste0(save_names[i], "_qcd.rds")))

save_bs_path_sample <- paste0(save_bs_path, foldername, "/", save_names[i], "/")
# vis: prep coords and count matrix
coords <- data.frame(x = spatialCoords(sce)[, 1], y = spatialCoords(sce)[, 2]); rownames(coords) <- colnames(sce)
counts_vis <- assay(sce, vis_rawmat)

ref <- Reference(chrom_mat, cell_types)
puck <- SpatialRNA(coords, counts_vis)
# RCTD --------------------------------------------------------------------
myRCTD <- create.RCTD(puck, ref, max_cores = 4)
myRCTD <- run.RCTD(myRCTD, doublet_mode = 'full')

norm_weights <- normalize_weights(myRCTD@results$weights)
RCTD_results <- data.frame(as(norm_weights, "matrix"))
write.csv(RCTD_results, paste0(save_bs_path_sample, save_names[i], "_spot_Level4_decon_RCTD.csv"))
results <- read.csv(paste0(save_bs_path_sample, save_names[i], "_spot_Level4_decon_RCTD.csv"), row.names = 1)


# Plot result -------------------------------------------------------------
library(magick)
image <- image_read(file.path(paste0(vis_path, "spatial"), "tissue_hires_image.png"))
info <- image_info(image)

scalef <- rjson::fromJSON(file = file.path(paste0(vis_path, "spatial"), "scalefactors_json.json"))

vis$pxl_col_in_fullres <- spatialCoords(vis)[, 1]
vis$pxl_row_in_fullres <- spatialCoords(vis)[, 2]

vis$pxl_col_in_hires <- vis$pxl_col_in_fullres * scalef$tissue_hires_scalef
vis$pxl_row_in_hires <- vis$pxl_row_in_fullres * scalef$tissue_hires_scalef

p0 <- ggplot(
  data.frame(x = 0, y = 0),
  aes(x, y)
) +
  geom_blank() +
  coord_fixed(
    expand = FALSE,
    xlim = c(min(vis$pxl_col_in_hires) - info$width * 0.05, max(vis$pxl_col_in_hires) + info$width * 0.05),
    ylim = c(min(vis$pxl_row_in_hires) - info$height * 0.05, max(vis$pxl_row_in_hires) + info$height * 0.05)
  ) +
  annotation_raster(
    image_flip(image),
    0,
    info$width,
    0,
    info$height
  ) +
  theme_bw()

library(randomcoloR)
set.seed(123)
colors <- distinctColorPalette(ncol(results))
radius <- scalef$spot_diameter_fullres * scalef$tissue_hires_scalef / 1.7

location <- data.frame(x = vis$pxl_col_in_hires,
                       y = vis$pxl_row_in_hires) 
rownames(location) <- colnames(vis)
colnames(location) <- c("x", "y")

overlap_barcodes <- intersect(rownames(results), rownames(location))
location <- location[overlap_barcodes, ]
results <- results[overlap_barcodes, ]
data = cbind(results, location)
ct.select = colnames(results)

library(scatterpie)
p <- p0 + 
  geom_scatterpie(aes(x = x, y = y, r = radius), 
                  data = data, cols = ct.select, color = NA) + 
  scale_fill_manual(values =  colors) +
  theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_rect(colour = "grey89", fill=NA, size=0.5),
        axis.text =element_blank(),
        axis.ticks =element_blank(),
        axis.title =element_blank(),
        legend.title=element_text(size = 16,face="bold"),
        legend.text=element_text(size = 15),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.key.size = unit(0.45, 'cm'),
        strip.text = element_text(size = 16,face="bold"),
        legend.position="right")+
  guides(fill=guide_legend(title="Cell Type", ncol = 1))


image_height <- (range(spatialCoords(vis)[, 2])[2] - range(spatialCoords(vis)[, 2])[1])
image_width <- (range(spatialCoords(vis)[, 1])[2] - range(spatialCoords(vis)[, 1])[1])

plot_width <- 20
plot_height <- plot_width * (image_height/image_width) - 2

fig_path <- here::here("./EuroBioc_results")
pdf(file = file.path(fig_path, "RCTD_deconHE.pdf"),
    width = plot_width,
    height = plot_height)
print(p)
dev.off()
