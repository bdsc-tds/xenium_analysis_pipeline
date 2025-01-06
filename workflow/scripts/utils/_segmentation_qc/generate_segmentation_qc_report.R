log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Get parameters for the report
pr <- list(
  segmentation_qc_path = snakemake@input[[1]],
  gene_panel = snakemake@params[["gene_panel"]]
)

rmarkdown::render(
  snakemake@params[["rmd_file"]], 
  params = pr,
  output_file = snakemake@output[[1]],
  clean = TRUE
)
