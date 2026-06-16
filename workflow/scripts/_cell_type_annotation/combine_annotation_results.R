log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

library(dplyr)
library(arrow)
library(stringr)
library(purrr)

input_dir <- snakemake@input[["root_dir"]]

# Find all labels.parquet files under input_dir
label_files <- list.files(
  path = input_dir,
  pattern = "^labels\\.parquet$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(label_files) == 0) {
  stop("No labels.parquet files found under: ", input_dir)
}


read_and_annotate <- function(path) {
  df <- read_parquet(path)
  
  # Extract tool/level/mode from path
  parts <- str_split(path, .Platform$file.sep)[[1]]
  n <- length(parts)
  
  if (n < 4) {
    stop("Unexpected path structure (too few components): ", path)
  }
  
  mode <- parts[n - 1]
  level <- parts[n - 2]
  tool  <- parts[n - 3]
  
  prefix <- paste(level, tool, mode, sep = "/")
  
  # Rename columns except cell_id
  colnames(df)[colnames(df) != "cell_id"] <- paste0(prefix, "/", colnames(df)[colnames(df) != "cell_id"])
  return(df)
}

merged <- label_files %>%
  map(read_and_annotate) %>%
  reduce(full_join, by = "cell_id")

write_parquet(merged, snakemake@output[["combined"]])
