"""
String constants for keys to be added to the configuration dictionary during processing to prepare for the execution of the workflow.
"""

# Child of the `config` dictionary (first level); for all wildcards.
WILDCARDS_NAME: str = "_wildcards"
# Child of `WILDCARDS_NAME` (second level); a list of values for diseases.
WILDCARDS_DISEASES_NAME: str = "_diseases"
# Child of `WILDCARDS_NAME` (second level); a list of values for gene panels.
WILDCARDS_GENE_PANELS_NAME: str = "_gene_panels"
# Child of `WILDCARDS_NAME` (second level); a list of values for donors.
WILDCARDS_DONORS_NAME: str = "_donors"
# Child of `WILDCARDS_NAME` (second level); a list of values for samples.
WILDCARDS_SAMPLES_NAME: str = "_samples"
# Child of `WILDCARDS_NAME` (second level); a list of values for segmentation methods.
WILDCARDS_SEGMENTATION_NAME: str = "_segmentation"
# Child of `WILDCARDS_NAME` (second level); a list of values for cell type annotation.
WILDCARDS_CELL_TYPE_ANNOTATION_NAME: str = "_cell_type_annotation"


# Child of `experiments` (second level); for the path to the experiments configuration file on the disk.
EXPERIMENTS_CONFIG_PATH_NAME: str = "_config_path"
# Child of `experiments` (second level); for the base path to the experiments on the disk.
EXPERIMENTS_BASE_PATH_NAME: str = "_base_path"
# Child of `experiments` (second level); for the name of or path to the gene panel files.
EXPERIMENTS_GENE_PANEL_FILES_NAME: str = "_gene_panel_file"
# Child of `experiments` (second level); for QC thresholds on gene panel level.
EXPERIMENTS_GENE_PANEL_QC_NAME: str = "_qc"
# Child of `experiments` (second level); colletions of samples on different levels.
EXPERIMENTS_COLLECTIONS_NAME: str = "_collections"
# Child of EXPERIMENTS_COLLECTIONS_NAME (third level); a collection of samples on disease level.
EXPERIMENTS_COLLECTIONS_DISEASES_NAME: str = "_diseases"
# Child of EXPERIMENTS_COLLECTIONS_NAME (third level); a collection of samples on gene panel level.
EXPERIMENTS_COLLECTIONS_GENE_PANELS_NAME: str = "_gene_panels"
# Child of EXPERIMENTS_COLLECTIONS_NAME (third level); a collection of samples on donors level.
EXPERIMENTS_COLLECTIONS_DONORS_NAME: str = "_donors"
# Child of `experiments` (second level); for cell type annotation on disease level.
EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME: str = "_cell_type_annotation"
# Child of EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME (third level); for paths to reference files for cell type annotation on disease level.
EXPERIMENTS_CELL_TYPE_ANNOTATION_REFERENCE_FILES_NAME: str = "paths"
# Child of EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME (third level); for resolution levels for cell type annotation on disease level.
EXPERIMENTS_CELL_TYPE_ANNOTATION_LEVELS_NAME: str = "levels"
# Child of EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME (third level); for the minimum number of cells required per cell type for cell type annotation on disease level.
EXPERIMENTS_CELL_TYPE_ANNOTATION_CELL_MIN_INSTANCES_NAME: str = "cell_min_instances"
