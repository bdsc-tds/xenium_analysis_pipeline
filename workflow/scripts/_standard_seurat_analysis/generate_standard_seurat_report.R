log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Get parameters for the report
pr <- list(
  xe_raw_path = snakemake@input[["raw"]],
  xe_path = snakemake@input[["preprocessed"]],
  default_assay = snakemake@params[["default_assay"]],
  segmentation_id = snakemake@params[["segmentation_id"]],
  sample_id = snakemake@params[["sample_id"]],
  normalisation_id = snakemake@params[["normalisation_id"]]
)

dir.create(
  file.path(snakemake@params[["intermediates_dir"]]),
  showWarnings = FALSE,
  recursive = TRUE,
  mode = "0766"
)
dir.create(
  file.path(snakemake@params[["knit_root_dir"]]),
  showWarnings = FALSE,
  recursive = TRUE,
  mode = "0766"
)

rmarkdown::render(
  snakemake@params[["rmd_file"]], 
  params = pr,
  output_file = snakemake@output[[1]],
  intermediates_dir = snakemake@params[["intermediates_dir"]],
  knit_root_dir = snakemake@params[["knit_root_dir"]],
  clean = TRUE
)

unlink(snakemake@params[["intermediates_dir"]], recursive = TRUE)
unlink(snakemake@params[["knit_root_dir"]], recursive = TRUE)
