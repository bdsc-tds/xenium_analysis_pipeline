######################################
#              Subrules              #
######################################

include: "_purification/purctd_full.smk" # depends on "_neighborhood_analysis/postprocess_rctd.smk" 
include: "_purification/purctd_balanced_by_spot_class.smk" # depends on "_purification/purctd_full.smk"
include: "_purification/purctd_balanced_by_score.smk" # depends on "_purification/purctd_full.smk" and _neighborhood_analysis/spatial_neighborhood_analysis.smk"
