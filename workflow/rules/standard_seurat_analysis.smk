import scripts._standard_seurat_analysis.standard_seurat_constants as sec


######################################
#              Subrules              #
######################################

include: "_standard_seurat_analysis/load_segmentation2seurat.smk"
include: "_standard_seurat_analysis/standard_qc.smk"
include: "_standard_seurat_analysis/standard_lognorm.smk"
include: "_standard_seurat_analysis/standard_sctransform_dimred_clust.smk"
include: "_standard_seurat_analysis/standard_seurat_report.smk"
