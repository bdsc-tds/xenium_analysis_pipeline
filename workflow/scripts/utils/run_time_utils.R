generate_xr_style_cell_names <- function(n, len = 8) {
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
