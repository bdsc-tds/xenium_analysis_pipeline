#######################################
#                Rules                #
#######################################

# Converts the same reference object used by `runSeurat2h5` (the .rds
# resolved by get_path2reference4reference_based_annotation) into a 10x-style
# HDF5 count matrix + a Parquet metadata table, so downstream Python steps
# (RCTD/scanpy/AnnData, etc.) can load the reference without going through
# R/rhdf5.
#
# NOTE: `script:` path is relative to *this* rule file's location — as
# written it assumes the same nesting depth as seurat.R
# ("../../../scripts/utils/seurat_to_h5_parquet.R"); adjust the leading
# "../" count if you place this .smk file elsewhere.

rule runConvertSeuratReferenceToH5:
    input:
        rds=f"{{reference_stem}}.rds"
    output:
        h5=protected(f"{{reference_stem}}_counts.h5"),
        parquet=protected(f"{{reference_stem}}_metadata.parquet")
    params:
        output_prefix=lambda wildcards, output: output.h5[: -len("_counts.h5")],
        assay="RNA",
        layer="counts",
        genome="unknown"
    resources:
        mem_mb=lambda wildcards, input, attempt: min(input.size_mb * attempt**2 * 50, 200000),
        retry_idx=lambda wildcards, attempt: attempt
    log:
        f"{{reference_stem}}_convertSeuratReferenceToH5.log"
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_reference_preprocessing/convert_seurat_2_h5.R"
