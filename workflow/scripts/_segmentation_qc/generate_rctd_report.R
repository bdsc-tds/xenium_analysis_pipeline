log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "output")
sink(log, type = "message")

# Get parameters for the report
pr <- list(
  file_path = snakemake@input[[1]],
  gene_panel_id = snakemake@params[["gene_panel_id"]]
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
