#######################################
#                Rules                #
#######################################

rule gatherGenePanelWiseSegmentedSamples:
    input:
        lambda wildcards: get_input2_or_params4gatherSegmentedSamples(
            wildcards,
            level=1,
        )
    output:
        protected(f'{config["output_path"]}/joint_scanpy_analysis/{{segmentation_id}}/{{gene_panel_id}}/gene_panel_wise_samples.json')
    params:
        lambda wildcards: get_input2_or_params4gatherSegmentedSamples(
            wildcards,
            level=1,
            for_input=False,
        )
    log:
        f'{config["output_path"]}/joint_scanpy_analysis/{{segmentation_id}}/{{gene_panel_id}}/logs/gatherGenePanelWiseSegmentedSamples.log'
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512
    script:
        "../../scripts/_joint_scanpy_analysis/gather_segmented_samples.py"

#TODO
rule runGenePanelWiseAnalysis:
    input:
        f'{config["output_path"]}/joint_scanpy_analysis/{{segmentation_id}}/{{gene_panel_id}}/gene_panel_wise_samples.json'
    output:
        protected(f'{config["output_path"]}/joint_scanpy_analysis/{{segmentation_id}}/{{gene_panel_id}}/umap.parquet')
    params:
        lambda wildcards: get_params4runAnalysis(
            wildcards,
            level=1,
        )
    log:
        f'{config["output_path"]}/joint_scanpy_analysis/{{segmentation_id}}/{{gene_panel_id}}/logs/runGenePanelWiseAnalysis.log'
    resources:
        mem_mb=lambda wildcards, attempt: min(
            get_mem_mb4runAnalysis(
                wildcards,
                attempt,
                level=1,
                multiplier=4,
            ),
            102400,
        ),
        slurm_partition=lambda wildcards: get_slurm_gpu_partition_name(
            wildcards,
        ) if _use_gpu() else "cpu",
        slurm_extra=get_slurm_extra
    container:
        config["containers"]["python_cuda"]
    shell:
        "mamba run -n general_cuda python3 workflow/scripts/_joint_scanpy_analysis/joint_embeding.py "
        "-i {input} "
        "-o {output} "
        "--n_comps {params[0][n_comps]} "
        "--n_neighbors {params[0][n_neighbors]} "
        "--metric {params[0][metric]} "
        "--min_dist {params[0][min_dist]} "
        "--min_counts {params[0][min_counts]} "
        "--min_features {params[0][min_features]} "
        "--max_counts {params[0][max_counts]} "
        "--max_features {params[0][max_features]} "
        "--min_cells {params[0][min_cells]} "
        "{params[0][use_gpu]} "
        "-l {log}"
