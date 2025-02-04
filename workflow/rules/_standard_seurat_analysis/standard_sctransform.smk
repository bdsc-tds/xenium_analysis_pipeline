#######################################
#                Rules                #
#######################################

"""
Comment:

- future_globals_maxSize: By default 10GB is used to deal with the [issue](https://github.com/satijalab/seurat/issues/1845) of "Global size exceeds maximum allowed size" when running Seurat (see the solution [here](https://satijalab.org/seurat/archive/v3.0/future_vignette.html)).
"""
rule runStandardScTransform:
    input:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/standard_qc/qced_seurat.rds'
    output:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/normalised_seurat.rds'
    params:
        future_globals_maxSize=lambda wildcards, resources: min(10**10 * resources[1], 10**11),
        default_assay=sec.SEURAT_DEFAULT_ASSAY,
        normalisation_id="sctransform"
    resources:
        mem_mb=lambda wildcards, input, attempt: max(input.size_mb * attempt * 50, 10240),
        retry_idx=lambda wildcards, attempt: attempt
    log:
        f'{config["output_path"]}/std_seurat_analysis/{{segmentation_id}}/{{sample_id}}/sctransform/logs/runStandardScTransform.log'
    container:
        config["containers"]["r"]
    script:
        "../../scripts/_standard_seurat_analysis/standard_sctransform.R"
