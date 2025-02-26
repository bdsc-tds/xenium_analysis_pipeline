#######################################
#              Functions              #
#######################################

def get_input2combineGenePanelWiseRCTDInfo(wildcards) -> list[str]:
    prefix: str = f'{config["output_path"]}/segmentation_qc'
    suffix: str = f'rctd_info.parquet'

    # Get all the annotations for the condition of the gene panel.
    annotations: list[str] = [
        i
        for i in get_dict_value(
            config,
            cc.WILDCARDS_NAME,
            cc.WILDCARDS_CELL_TYPE_ANNOTATION_NAME,
            extract_layers_from_experiments(
                wildcards.gene_panel_id,
                0,
            )[0],
        )
        if re.match(
            r".+/rctd_.+",
            i,
            flags=re.IGNORECASE,
        ) is not None
    ]

    # Get all the samples for the gene panel.
    samples: list[str] = get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_COLLECTIONS_NAME,
        cc.EXPERIMENTS_COLLECTIONS_GENE_PANELS_NAME,
        wildcards.gene_panel_id,
    )

    return [
        os.path.join(
            prefix,
            i,
            sample,
            j,
            k,
            suffix,
        )
        for i in SEGMENTATION_ID
        for sample in samples
        for j in NORMALISATION_ID
        for k in annotations
    ]


#######################################
#                Rules                #
#######################################

rule combineGenePanelWiseRCTDInfo:
    input:
        get_input2combineGenePanelWiseRCTDInfo
    output:
        protected(f'{config["output_path"]}/segmentation_qc/{{gene_panel_id}}/rctd_info.parquet')
    log:
        f'{config["output_path"]}/segmentation_qc/{{gene_panel_id}}/logs/combineGenePanelWiseRCTDInfo.log'
    conda:
        "../../envs/pyarrow.yml"
    resources:
        mem_mb=lambda wildcards, input, attempt: min(input.size_mb * attempt * 2000, 20480)
    script:
        "../../scripts/_segmentation_qc/combine_rctd_info.py"
