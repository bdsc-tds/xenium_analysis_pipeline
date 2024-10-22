######################################
#              Subrules              #
######################################

include: '_standard_seurat_analysis/convert_segmentation_to_seurat.smk'
include: '_standard_seurat_analysis/standard_QC_xenium.smk'
include: '_standard_seurat_analysis/standard_lognorm_xenium.smk'
include: '_standard_seurat_analysis/standard_SCT_dimred_clust_xenium.smk'
include: '_standard_seurat_analysis/generate_report.smk'

