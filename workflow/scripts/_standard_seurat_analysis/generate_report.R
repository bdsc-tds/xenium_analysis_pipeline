log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output") # @Senbai, is it fine to change output stream when generating html report? 
sink(log, type = "message")

# Get parameters for the report
pr <- list(
  xe_raw_path    <- snakemake@input[["xe_raw"]],
  xe_path        <- snakemake@input[["xe"]],
  default_assay  <- snakemake@params[["default_assay"]],
)

rmarkdown::render(
  "report.Rmd", 
  params = pr,
  output_file = snakemake@output[["report"]],
  clean = TRUE
)