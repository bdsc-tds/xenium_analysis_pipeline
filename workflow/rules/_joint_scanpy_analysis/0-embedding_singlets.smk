# cfg paths
xenium_processed_data_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])

# params from pipeline config
normalisations = ['lognorm','sctransform']
layers = ['data','scale_data']
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# Params
n_comps = config['umap_n_comps']
n_neighbors = config['umap_n_neighbors']
min_dist = config['umap_min_dist']
metric = config['umap_metric']

# to get singlet info
references = ['matched_reference_combo','external_reference']
methods = ['rctd_class_aware']
levels = ['Level2.1',]#['Level1','Level2','Level3','Level4',] # condition and sample as color to plot added here in addition to levels


out_files_panel = []

for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem == 'proseg_mode':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation in normalisations: 
                for layer in layers:
                    for reference in references:
                        for method in methods:
                            for level in levels:
                                if level == 'Level2.1' and reference == 'external_reference':
                                    continue

                                k = (segmentation.stem,condition.stem,panel.stem,normalisation,reference,method,level)
                                name = '/'.join(k)
                                rule_name = '/'.join(k+(layer,))

                                out_file = results_dir / f'embed_panel_singlets/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                                out_files_panel.append(out_file)

                                rule:
                                    name: f'embed_panel_singlets/{rule_name}'
                                    input:
                                        panel=panel,
                                    output:
                                        out_file=out_file,
                                    params:
                                        xenium_processed_data_dir=xenium_processed_data_dir,
                                        normalisation=normalisation,
                                        layer=layer,
                                        n_comps=n_comps,
                                        n_neighbors=n_neighbors,
                                        metric=metric,
                                        min_dist=min_dist,
                                        min_counts=min_counts,
                                        min_features=min_features,
                                        max_counts=max_counts,
                                        max_features=max_features,
                                        min_cells=min_cells,
                                        cell_type_annotation_dir=cell_type_annotation_dir,
                                        reference=reference,
                                        method=method,
                                        level=level,
                                    threads: 1
                                    resources:
                                        mem='100GB' if panel.stem == '5k' else '50GB',
                                        # runtime='30m' if panel.stem == '5k' else '20m',
                                        runtime='8h' if panel.stem == '5k' else '3h',
                                        # slurm_partition = "gpu",
                                        # slurm_extra = '--gres=gpu:1',
                                    conda:
                                        "spatial"
                                    shell:
                                        """
                                        mkdir -p "$(dirname {output.out_file})"

                                        python workflow/scripts/xenium/embed_panel_singlets.py \
                                            --panel {input.panel} \
                                            --out_file {output.out_file} \
                                            --xenium_processed_data_dir {params.xenium_processed_data_dir} \
                                            --normalisation {params.normalisation} \
                                            --layer {params.layer} \
                                            --n_comps {params.n_comps} \
                                            --n_neighbors {params.n_neighbors} \
                                            --metric {params.metric} \
                                            --min_dist {params.min_dist} \
                                            --min_counts {params.min_counts} \
                                            --min_features {params.min_features} \
                                            --max_counts {params.max_counts} \
                                            --max_features {params.max_features} \
                                            --min_cells {params.min_cells} \
                                            --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                                            --reference {params.reference} \
                                            --method {params.method} \
                                            --level {params.level} \

                                        echo "DONE"
                                        """


rule embed_panel_singlets_all:
    input:
        out_files_panel
