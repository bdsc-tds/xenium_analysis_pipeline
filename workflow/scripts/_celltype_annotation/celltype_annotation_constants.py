"""
String constants for reference (scRNA-seq) and query (Xenium) Seurat objects.
"""

REF_SEURAT_DEFAULT_ASSAY: str = "RNA"
XE_SEURAT_DEFAULT_ASSAY: str = "Xenium"
#XE_SEURAT_SPATIAL_DIM_NAME: str = "spatial" # Not used 

# The following parameters may vary depending on the size of the Xenium gene panel
# Filtering xenium cells for annotation (mb redundant, but if QC differs, its better to filter low count cells from annotation results)
XE_MIN_UMI: int = 10 # min total UMIs
XE_MIN_counts: int = 10 # min counts (ie., UMIs), when restricted to common genes

# Filtering reference cells 
REF_MIN_UMI: int = 10 # low, given the small number of counts in scRNA-seq data when restricted to ~300 plex xenium panel. may be usedful to be dependent on `gene_panel`, larger the panel -> more common genes between the ref and query data -> the more UMI expected to be present in the reference
REF_MAX_UMI: int = 2000 # low, given the small number of counts in scRNA-seq data when restricted to ~300 plex xenium panel.

# RCTD-specific constants
UMI_min_sigma: int = 100 # min UMI in xenium cells to consider cell for sigma extimation (sigma of  Log-Poisson distribution) #TODO: replace with one from `config -> ... rctd -> mode -> other_options`
CELL_MIN_INSTANCE: int = 25 # minimum number of cells required per cell type. Can be decreased for smaller references

# SingleR-specific params 
singleR_genes: str = "de" # use differentially expressed genes 
singleR_de_method: str = "t" # test to use for DEA
singleR_aggr_ref: bool = True 
singleR_aggr_args: dict = {"rank": 50, "power": 0.7}

# XGBoost-specific params
xgb_nrounds: int = 1000 
xgb_eta: float = 0.3 

# Seurat_Transfer-specific params
dims: list = list(range(1, 51))
