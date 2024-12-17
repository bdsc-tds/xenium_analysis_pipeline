"""
String constants for reference (scRNA-seq) and query (Xenium) Seurat objects.
"""

REF_SEURAT_DEFAULT_ASSAY: str = "RNA"
XE_SEURAT_DEFAULT_ASSAY: str = "Xenium"
# XE_SEURAT_SPATIAL_DIM_NAME: str = "spatial" # Not used

# The following parameters may vary depending on the size of the Xenium gene panel
# Filtering xenium cells for annotation (mb redundant, but if QC differs, its better to filter low count cells from annotation results)
XE_MIN_UMI: int = 10  # min total UMIs
XE_MIN_counts: int = 10  # min counts (ie., UMIs), when restricted to common genes

# Filtering reference cells
REF_MIN_UMI: int = (
    10  # low, given the small number of counts in scRNA-seq data when restricted to ~300 plex xenium panel. may be usedful to be dependent on `gene_panel`, larger the panel -> more common genes between the ref and query data -> the more UMI expected to be present in the reference
)

REF_MAX_UMI: int = (
    2000  # low, given the small number of counts in scRNA-seq data when restricted to ~300 plex xenium panel.
)
