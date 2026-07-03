######################################
#             Functions              #
######################################


def get_path2reference_by_reference_path(wildcards) -> str:
    """
    Resolve the reference_path wildcard (e.g. "CRC/matched_reference") back
    to its source .rds file, using the table built in the wildcards file.
    """
    return REFERENCE_PATH_TO_RDS[wildcards.reference_path]
  

######################################
#              Subrules              #
######################################

include: '_reference_preprocessing/convert_seurat_2_h5.smk'
