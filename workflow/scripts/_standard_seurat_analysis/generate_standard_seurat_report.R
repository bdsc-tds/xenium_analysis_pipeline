log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Get parameters for the report
pr <- list(
  xe_raw_path = snakemake@input[["raw"]],
  xe_path = snakemake@input[["preprocessed"]],
  default_assay = snakemake@params[["default_assay"]],
  segmentation_id = snakemake@params[["segmentation_id"]],
  sample_id = snakemake@params[["sample_id"]]
)

rmarkdown::render(
  snakemake@params[["rmd_file"]], 
  params = pr,
  output_file = snakemake@output[[1]],
  clean = TRUE
)