library(Matrix)
library(SeuratObject)
library(Seurat)
library(dplyr)

generate_xr_style_cell_names <- function(n, len = 8, seed = 42) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    unique_strings <- character(0)
    while (length(unique_strings) < n) {
        new_string <- paste0(
            paste0(
                sample(letters, len, replace = TRUE),
                collapse = ""
            ),
            "-1"
        )

        if (!new_string %in% unique_strings) {
            unique_strings <- c(unique_strings, new_string)
        }
    }

    return(unique_strings)
}

# Change 1-based index to 0-based 
quick_generate_cell_names <- function(n, prefix = "") {
    return(paste0(prefix, as.integer(seq(n) - 1)))
}

# Assume `new_counts` is cell-by-gene.
# Assume `cell_coords` is, if provided, a data frame with two columns for x and y coordinates.
replace_counts_in_seurat <- function(xe, new_counts, cell_coords = NULL, cell_id_prefix = "") {
    generate_cell_ids <- sum(!grepl("\\w{8}-1", rownames(new_counts))) > 0

    if (generate_cell_ids) {
        # if (nrow(new_counts) < 1e5) {
        #     cell_ids <- generate_xr_style_cell_names(nrow(new_counts))
        # } else {
            cell_ids <- quick_generate_cell_names(nrow(new_counts), cell_id_prefix)
        # }
    } else {
        cell_ids <- rownames(new_counts)
    }

    if (is.null(cell_coords)) {
        if (generate_cell_ids) {
            stop("Error! Missing cell coordinates, or the count matrix must have rownames from the provided Seurat object.")
        }

        cell_coords <- data.frame(
            cell = cell_ids
        ) %>%
            inner_join(
                GetTissueCoordinates(xe),
                by = "cell"
            ) %>%
            select(
                -cell
            )
        
        if (nrow(cell_coords) != length(cell_ids)) {
            stop("Error! The counts matrix has cells not existing in the passed Seurat object.")
        }
    } else {
        if (!is.data.frame(cell_coords)) {
            stop("Error! Cell coordinates must be passed as a data frame.")
        }

        if (nrow(new_counts) != nrow(cell_coords)) {
            stop("Error! Different cell numbers detected in the count matrix and cell coordinates.")
        }

        if (ncol(cell_coords) != 2) {
            stop("Error! The data frame for cell coordinates should only have two columns for x and y coordinates.")
        }
    }

    gene_ids <- colnames(new_counts)

    new_counts <- t(as.matrix(new_counts))

    # Create a new Seurat object for expected counts
    new_xe <- CreateSeuratObject(
        counts = new_counts,
        assay = "Xenium"
    )
    colnames(new_xe) <- cell_ids
    rownames(new_xe) <- gene_ids

    coords <- CreateCentroids(
        coords = cell_coords,
        nsides = xe@images$fov@boundaries$centroids@nsides,
        radius = xe@images$fov@boundaries$centroids@radius,
        theta = xe@images$fov@boundaries$centroids@theta
    )
    coords <- RenameCells(coords, cell_ids)

    new_xe@images <- list(
        fov = CreateFOV(
        coords,
        molecules = NULL,
        assay = xe@images$fov@assay,
        key = xe@images$fov@key
        )
    )

    missing_assays <- Assays(xe)[!Assays(xe) %in% Assays(new_xe)]
    names(missing_assays) <- missing_assays

    new_assays <- sapply(
        missing_assays,
        function(x) {
            assay <- GetAssayData(xe, assay = x)

            CreateAssayObject(
                counts = Matrix(
                    0,
                    ncol = ncol(new_xe),
                    nrow = nrow(assay),
                    dimnames = list(
                        rownames(assay),
                        cell_ids
                    )
                )
            )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    for (i in names(new_assays)) {
        new_xe[[i]] <- new_assays[[i]]
    }

    return(new_xe)
}
