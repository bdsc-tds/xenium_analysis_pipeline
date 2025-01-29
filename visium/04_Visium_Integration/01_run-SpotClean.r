library(Seurat)
library(SpotClean)
library(SummarizedExperiment)
library(cowplot)

keepHighGene =function(count_mat, top_high=5000,
                         mean_cutoff=1, return_matrix=FALSE,
                         verbose=TRUE){

    mean_exp <- rowMeans(count_mat)
    # keep at most this number of highly expressed genes.
    #top_genes <- rank(-mean_exp)<=top_high in SpotClean code but not convinced
    top_genes <- order(mean_exp,decreasing=TRUE)[1:min(top_high, length(mean_exp))]## add a condition to take ngenes in count_mat if nrow < top_high
    count_mat <- count_mat[top_genes,,drop=FALSE]
    mean_exp <- mean_exp[top_genes]

    high_exp_genes <- mean_exp>=mean_cutoff

    S_vf <- NormalizeData(CreateSeuratObject(count_mat), verbose = FALSE)
    S_vf <- FindVariableFeatures(S_vf,
                                 selection.method = "mvp", verbose = FALSE)
    
    # accommodate changes in Seurat v5 objects
    if(as.integer(gsub("\\<(\\d+)\\.\\d+\\.\\d+", "\\1", S_vf@version))>=5){
        high_variable_genes <- S_vf@assays$RNA@meta.data$vf_mvp_data_variable
    }else{
        high_variable_genes <- S_vf@assays$RNA@meta.features$mvp.variable
    }
    
    gene_tokeep <- high_variable_genes | high_exp_genes

    if(verbose){
        message("Kept ",sum(gene_tokeep),
                " highly expressed or highly variable genes.")
    }

    if(return_matrix){
        return(count_mat[gene_tokeep,, drop=FALSE])
    }else{
        return(names(which(gene_tokeep)))
    }

}

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

datapath.spcl <- paste0("/work/PRTNR/CHUV/DIR/rgottar1/spatial/Owkin_Pilot_Data/Visium_integration_rep_owkin/Seurat5_SpCl1.4.1_final/", disease, "/spotclean/")


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

  
  # Load Visium data
  sample_raw = read10xRaw(file.path(data_path,"raw_feature_bc_matrix"))
  sample_slide_info = read10xSlide(paste0(data_path, "spatial/tissue_positions.csv"), 
                                   paste0(data_path, "spatial/tissue_lowres_image.png"), 
                                   paste0(data_path, "spatial/scalefactors_json.json"))
  
  ## sample_raw and sample_slide_info must have the same barcodes
  sample_raw = sample_raw[,intersect(sample_slide_info$slide$barcode,colnames(sample_raw))]
  sample_slide_info$slide = sample_slide_info$slide[sample_slide_info$slide$barcode %in% intersect(sample_slide_info$slide$barcode,colnames(sample_raw)),]
  
  ## genes not in the filtered probe set have to be removed from the sample_raw (probe_set csv file, included == TRUE), we keep the rownames that are common with the filtered matrix
  sample_filtered = read10xRaw(paste0(data_path, "filtered_feature_bc_matrix")) 
  sample_raw = sample_raw[rownames(sample_filtered), ]
  
  # Visualize raw data
  sample_obj <- createSlide(count_mat = sample_raw, 
                            slide_info = sample_slide_info)
  
  metadata(sample_obj)$slide$total_counts <- Matrix::colSums(sample_raw)
  
  tissue.obj=sample_obj[, which(metadata(sample_obj)$slide$tissue==1)]
  tissue.obj@metadata$slide = sample_obj@metadata$slide[which(metadata(sample_obj)$slide$tissue==1), ]
  # Decontaminate raw data
  # the keepHighgene is implementing in Spotclean but bugging for now... so I make it in 2 steps
  gene_keep <- keepHighGene(assays(sample_obj)$raw[,which(metadata(sample_obj)$slide$tissue==1),  drop=FALSE], verbose=TRUE)
  decont_obj <- spotclean(sample_obj,gene_keep=gene_keep)
  
  # (Optionally) Transform to Seurat object for downstream analyses
  seurat_obj <- convertToSeurat(decont_obj, image_dir = paste0(data_path, "spatial/"))
  saveRDS(seurat_obj, file=paste0(datapath.spcl, "/SpotClean_", samplename, ".rds"))
  ## Create the matrix of corrected counts in csv format
}


for(n in 1:5){
  samplename <-sample_name[n]
  print(dim(readRDS(paste0(datapath.spcl, "SpotClean_", samplename, ".rds"))))
}



