#' Extract Parameters from a Reference-Based annotation_id
#'
#' This function parses a reference-based annotation ID string and extracts parameters based on the expected components.
#' It checks that the annotation ID is of type "reference_based" and that the ID has the correct number of components to match 
#' the provided parameter names.
#'
#' @param annotation_id A character string representing the annotation ID, where components are separated by "/".
#'                      The ID should begin with "reference_based" and contain additional parts that represent specific parameters.
#' @param parameters_names A character vector of parameter names to assign to each component in `annotation_id`. 
#'                         Defaults to `c("annotation_approach", "reference_type", "annotation_method", 
#'                         "annotation_level", "annotation_mode")`.
#' 
#' @return A named character vector where each component from `annotation_id` is assigned to a corresponding parameter 
#'         from `parameters_names`.
#'
#' @examples
#' # Example usage:
#' annotation_id <- "reference_based/type1/method2/level3/mode4"
#' get_referencebased_annotation_parameters(annotation_id)
#'
#' @throws If `annotation_id` does not start with "reference_based" or if the number of components in `annotation_id` does not
#'         match the length of `parameters_names`, an error is thrown.

get_referencebased_annotation_parameters <- function(
    annotation_id,
    parameters_names = c("annotation_approach", "reference_type", "annotation_method", "annotation_level", "annotation_mode")
) {
  # Split the annotation_id string into components
  res <- strsplit(annotation_id, "/")[[1]]
  
  # Check if the first component is "reference_based"
  if (res[1] != "reference_based") {
    stop(sprintf("`reference_based` annotation is expected, but %s provided", res[1]))
  }
  
  # Check if the number of components matches the length of parameters_names
  if (length(res) != length(parameters_names)) {
    stop("Non-matching number of annotation parameters in `annotation_id` and `parameters_names`")
  }
  
  # Assign names to each component based on parameters_names
  names(res) <- parameters_names
  
  return(res)
}


#' Generate a reference object by filtering and subsetting cells
#'
#' This function creates a reference object by filtering and subsetting cells in 
#' a single-cell dataset (`chrom`) according to specified parameters. The function 
#' allows for matched or external reference generation, annotating at various levels, 
#' and filtering cells based on gene expression metrics.
#'
#' @param chrom A `Seurat` object containing the input data (e.g., single-cell RNA-seq data).
#' @param query_features Character vector. Vector of query features (\code{rownames(xe)})
#' @param donor_id Character. Query donor ID
#' @param reference_type Character. Specify "matched" for donor-specific reference or "external" for a general reference.
#' @param annotation_level Character. Level of cell annotation used for subsetting (e.g., "Level3").
#' @param ref_assay Character. Assay to use in the reference object (e.g., "RNA").
#' @param MIN_UMI_ref Numeric. Minimum UMI count for filtering low-expression cells.
#' @param MAX_UMI_ref Numeric. Maximum UMI count for filtering high-expression cells (potential doublets).
#' @param CELL_MIN_INSTANCE Numeric. Minimum number of cells per cell type to retain the type.
#'
#' @return A filtered `Seurat` object with cells meeting specified criteria.
#'
generate_reference_obj <- function(
    chrom,
    query_features,
    donor_id,
    # reference_type = "matched_reference", # deprecated 
    annotation_level = "Level3",
    ref_assay = "RNA",
    MIN_UMI_ref = 10,
    MAX_UMI_ref = 2000,
    CELL_MIN_INSTANCE = 25
) {
  
  # Check if donor is in the reference 
  if ("donor" %in% colnames(chrom@meta.data)) {
    ref_donors <- unique(chrom@meta.data$donor)
    DONOR_IN_REF <- donor_id %in% ref_donors
  } else {
    DONOR_IN_REF <- FALSE
  }
  
  # Filter non-tumor cells or donor-specific cells based on reference type
  if (DONOR_IN_REF) {
    print(donor_id)
    orig_cell_type <- unique(chrom[[annotation_level]] %>% pull())
    orig_ncells <- ncol(chrom)
    
    chrom$isDonorSpecific <- grepl(paste(ref_donors, collapse = "|"), chrom[[annotation_level]] %>% pull())
    chrom <- subset(chrom, subset = ((isDonorSpecific == FALSE) | (donor == donor_id)))
    
    message(sprintf("For the reference dataset, %d cells were removed from the original dataset.", orig_ncells - ncol(chrom)))
    message(sprintf("The following cell types were removed: %s", paste(setdiff(orig_cell_type, unique(chrom[[annotation_level]] %>% pull())), collapse = ", ")))
  }
  
  all_cell_types <- as.factor(chrom@meta.data %>% pull(annotation_level)) %>% levels()
  
  # Subset reference data to used reference and panel genes
  chrom <- DietSeurat(chrom, assays = ref_assay, dimreducs = NULL, graphs = NULL)
  #chrom <- subset(chrom, features = rownames(chrom)[rownames(chrom) %in% query_features])
  
  features <- rownames(chrom)[rownames(chrom) %in% query_features]
  # Filter cells based on UMI counts to exclude low and high count cells
  chrom$nCount <- colSums(GetAssayData(chrom, assay = ref_assay, layer = "counts")[features,]) #chrom@meta.data[[paste0("nCount_", ref_assay)]]
  chrom <- subset(chrom, subset = (nCount > MIN_UMI_ref & nCount < MAX_UMI_ref))
  
  # Remove cell types with fewer than CELL_MIN_INSTANCE cells
  ref_labels <- as.factor(chrom@meta.data %>% pull(annotation_level))
  names(ref_labels) <- colnames(chrom)
  keep_cell_types <- names(which(summary(ref_labels) >= CELL_MIN_INSTANCE))
  ref_labels <- ref_labels[ref_labels %in% keep_cell_types] %>% droplevels()
  chrom <- subset(chrom, cells = names(ref_labels))
  
  if(length(all_cell_types) != length(keep_cell_types)){
    message(sprintf("Cell types %s were filtered out from the reference due to low abundance and/or low profile coverage.", 
                    paste(setdiff(all_cell_types, keep_cell_types), collapse = ",")))
  }
  
  return(chrom)
}



#' Find top positive markers and return binary celltype x gene matrix
#'
#' @param obj Either a Seurat object or a gene expression matrix (genes x cells)
#' @param celltypes A vector of cell type labels (length = ncol(obj) if matrix)
#' @param top_n Maximum number of markers to select per cell type
#' @param min_padj Adjusted p-value cutoff
#'
#' @return A binary matrix (celltypes x genes)
generate_marker_binary_matrix <- function(obj, celltypes, top_n = 20, min_padj = 0.05) {
  
  # Ensure input is Seurat
  if (inherits(obj, "Seurat")) {
    seu <- obj
  } else if (is.matrix(obj)) {
    if (is.null(colnames(obj))) colnames(obj) <- paste0("cell", seq_len(ncol(obj)))
    seu <- Seurat::CreateSeuratObject(counts = obj)
  } else {
    stop("Input must be a Seurat object or expression matrix")
  }
  
  # Add cell type metadata
  if (length(celltypes) != ncol(seu)) {
    stop("Length of celltypes vector must match number of cells")
  }
  seu$cluster <- factor(celltypes)
  
  Seurat::Idents(seu) <- "cluster"
  # DEA per cell type
  all_markers <- Seurat::FindAllMarkers(
    object = seu, 
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  
  # Filter by adjusted p-value
  all_markers <- all_markers[all_markers$p_val_adj < min_padj, ]
  
  # For each cluster, take top N genes
  top_markers <- all_markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  
  # Build binary matrix
  celltypes_unique <- unique(top_markers$cluster)
  genes_unique <- unique(top_markers$gene)
  bin_mat <- matrix(0, nrow = length(celltypes_unique), ncol = length(genes_unique),
                    dimnames = list(celltypes_unique, genes_unique))
  
  for (i in seq_len(nrow(top_markers))) {
    bin_mat[top_markers$cluster[i], top_markers$gene[i]] <- 1
  }
  
  return(bin_mat)
}



#' Generate a `class_df` Data Frame for Annotation
#'
#' This function creates a `class_df` data frame by selecting and processing 
#' specified annotation levels from the metadata of a `Seurat` object. If the 
#' `class_level` annotation is valid and present, the function constructs a unique 
#' mapping between `annotation_level` and `class_level`. The resulting `class_df` 
#' is intended for use in cell classification tasks.
#'
#' @param chrom A `Seurat` object containing metadata with cell annotation information.
#' @param annotation_level Character. The primary annotation level to use as row names in `class_df`.
#' @param class_level Character. The secondary annotation level to be mapped to `annotation_level`.
#'
#' @return A data frame (`class_df`) with `annotation_level` as row names and `class` as a column.
#' Returns `NULL` if `class_level` is not specified or not found in `chrom@meta.data`.
#'
#' @examples
#' # Example usage:
#' # class_df <- generate_class_df(chrom = my_seurat_obj, annotation_level = "Level3", class_level = "Level2")
#'
generate_class_df <- function(
    chrom,
    annotation_level,
    class_level
) {
  # Initialize class_df as NULL
  class_df <- NULL
  
  if (!is.null(class_level) && class_level == annotation_level) {
    class_level <- NULL
  }
  
  # Check if class_level is specified and present in chrom metadata
  if (!is.null(class_level)) {
    if (class_level %in% names(chrom@meta.data)) {
      # Generate class_df with distinct annotation-class mappings
      class_df <- chrom@meta.data %>%
        select(all_of(c(annotation_level, class_level))) %>%
        group_by(across(all_of(c(annotation_level, class_level)))) %>%
        summarise(n = n(), .groups = "drop") %>%
        group_by(across(all_of(annotation_level))) %>%
        slice_max(n, with_ties = FALSE) %>% # Keep the most frequent class for each annotation_level
        ungroup() %>%
        rename(class = all_of(class_level)) %>%
        select(-n) %>% # Remove the count column
        as.data.frame() # Convert to base R data.frame
      
      # Set annotation_level as row names
      rownames(class_df) <- class_df[[annotation_level]]
      class_df <- class_df[setdiff(names(class_df), annotation_level)]
    } else {
      warning(sprintf("`class_level` '%s' does not exist in the reference meta.data; RCTD will be run with `class_df = NULL`", class_level))
    }
  } else {
    warning("`class_level` was not provided; RCTD will be run with `class_df = NULL`")
  }
  
  return(class_df)
}


#'Convert a List of Matrices to a Long-Format Data Frame
#'
#' This function takes a list of matrices (mat_score) from rctd output and converts each matrix into a long-format data frame.
#' The resulting data frames are combined into a single data frame with an additional column
#' identifying the index of the source matrix in the input list.
#'
#' @param score_mat rctd@results$score_mat
#'
#' @return A data frame with the following columns:
#' \describe{
#'   \item{row}{Row indices of the matrix entries.}
#'   \item{column}{Column indices of the matrix entries.}
#'   \item{value}{Values from the matrices.}
#'   \item{cell_index}{The index of the matrix in the original list.}
#' }
#'
#' @examples
#' # Example usage:
#' score_mat <- convert_score_mat_to_long_df(rctd@results$score_mat)
#' print(score_mat)
convert_score_mat_to_long_df <- function(
    score_mat
){
  # Convert each matrix into a long-format data frame
  res <- lapply(seq_along(score_mat), function(idx) {
    mat <- score_mat[[idx]]
    df <- as.data.frame(as.table(as.matrix(mat)))
    df$cell_index <- idx  # Add the index as a new column
    return(df)
  })
  res <- data.table::rbindlist(res)
  colnames(res) <- c("row", "column", "value", "cell_index")

  return(res)
}


#' Convert Singlet Scores to a Long Data Frame
#'
#' This function converts a list of singlet scores from rctdf results (`rctd@results$singlet_scores`) into a long-format data frame, with each
#' score linked to its corresponding cell index and cell type.
#'
#' @param singlet_list a `rctd@results$singlet_scores`
#'
#' @return A tibble with the following columns:
#' \describe{
#'   \item{cell_index}{The index of the cell in the input list.}
#'   \item{cell_type}{The name of the cell type corresponding to each score.}
#'   \item{score}{The score value.}
#' }
#'
#' @examples
#' # Example usage
#' singlet_scores <- convert_singlet_scores_to_long_df(rctd@results$singlet_scores)
#' print(singlet_scores)
#'
convert_singlet_scores_to_long_df <- function(
    singlet_scores
){
  singlet_unlist <- unlist(singlet_scores, use.names = TRUE)
  singlet_tibble <- tibble(
    cell_index = lapply(seq_along(singlet_scores), function(idx) rep(idx, length(singlet_scores[[idx]]))) %>% unlist(),
    cell_type = names(singlet_unlist),          # Extract names
    score = singlet_unlist  # Extract values
  )
  return(singlet_tibble)
}
