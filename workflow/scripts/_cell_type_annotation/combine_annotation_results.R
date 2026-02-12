#!/usr/bin/env Rscript

library(optparse)
library(dplyr)
library(arrow)
library(stringr)
library(purrr)

# Argument parser
option_list <- list(
  make_option(c("--input_dir"), type="character", help="Input root directory"),
  make_option(c("--output_file"), type="character", help="Output file path")
)
opt <- parse_args(OptionParser(option_list=option_list))

# Find all labels.parquet files under input_dir
label_files <- list.files(
  path = opt$input_dir,
  pattern = "labels\\.parquet$",
  full.names = TRUE,
  recursive = TRUE
)

if (length(label_files) == 0) {
  stop("No labels.parquet files found.")
}

read_and_annotate <- function(path) {
  df <- read_parquet(path)
  
  # Extract tool/level/mode from path
  parts <- str_split(path, .Platform$file.sep)[[1]]
  n <- length(parts)
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

write_parquet(merged, opt$output_file)
