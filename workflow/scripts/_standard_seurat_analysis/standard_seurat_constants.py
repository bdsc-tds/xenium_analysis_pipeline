"""
String constants for Seurat objects.
"""

XENIUM_CONTROL_GENE_PAT: str = "^BLANK_\\d+|^UnassignedCodeword-\\d+|^NegControl.+"

SEURAT_SPATIAL_DIM_NAME: str = "spatial"

SEURAT_DEFAULT_ASSAY: str = "Xenium"
SEURAT_ALT_ASSAY: str = "SCT"

SEURAT_DEFAULT_LAYER: str = "counts"

SEURAT_DATA_LAYER: str = "data"
SEURAT_SCALE_DATA_LAYER: str = "scale.data"
