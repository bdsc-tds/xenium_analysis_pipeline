import scripts._celltype_annotation.celltype_annotation_constants.py as sac
include: "../../../scripts/_celltype_annotation/annotation_utils.py" # @Senbai, can I import funcrions that will be used by all listed below .smk like this?


######################################
#              Subrules              #
######################################

include: "_celltype_annotation/_reference_based/rctd.smk"
include: "_celltype_annotation/_reference_based/singleR.smk"
include: "_celltype_annotation/_reference_based/xgboost.smk"
include: "_celltype_annotation/_reference_based/seurat_transfer.smk"
#include: "_celltype_annotation/_reference_based/tangram.smk" # placeholder
