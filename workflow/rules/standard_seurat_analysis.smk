import scripts._standard_seurat_analysis.standard_seurat_constants as sec


######################################
#              Subrules              #
######################################

include: "_standard_seurat_analysis/load_segmentation2seurat.smk"
include: "_standard_seurat_analysis/standard_qc.smk"
include: '_standard_seurat_analysis/standard_lognorm_xenium.smk'
include: '_standard_seurat_analysis/standard_SCT_dimred_clust_xenium.smk'
include: '_standard_seurat_analysis/generate_report.smk'

