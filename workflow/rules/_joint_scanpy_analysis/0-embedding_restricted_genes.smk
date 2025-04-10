import pandas as pd

# cfg paths
xenium_processed_data_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])

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

# genes = pd.read_csv(config['markers_dir']+'Xenium_hLung_v1_metadata.csv')['Gene'].tolist()
genes = pd.read_csv(config['markers_dir']+'Xenium_NSCLC_5k_lung_chromium_common_genes.csv')['gene'].tolist()

out_files_panel = []

for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
    if segmentation.stem == 'proseg_mode':
        continue
    for condition in (conditions := segmentation.iterdir()): 
        for panel in (panels := condition.iterdir()):
            for normalisation in normalisations: 
                for layer in layers:
                    k = (segmentation.stem,condition.stem,panel.stem,normalisation)
                    name = '/'.join(k)
                    rule_name = '/'.join(k+(layer,))

                    out_file = results_dir / f'embed_panel_restricted_genes/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
                    out_files_panel.append(out_file)

                    rule:
                        name: f'embed_panel_restricted_genes/{rule_name}'
                        input:
                            panel=panel,
                        output:
                            out_file=out_file,
                        params:
                            xenium_processed_data_dir = xenium_processed_data_dir,
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
                            genes=genes,
                        threads: 1
                        resources:
                            mem='100GB' if panel.stem == '5k' else '50GB',
                            # runtime='30m' if panel.stem == '5k' else '20m',
                            runtime='12h',
                            # slurm_partition = "gpu",
                            # slurm_extra = '--gres=gpu:1',
                        conda:
                            "spatial"
                        shell:
                            """
                            mkdir -p "$(dirname {output.out_file})"

                            python -u workflow/scripts/xenium/embed_panel.py \
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
                                --genes {params.genes}
                                
                            echo "DONE"
                            """


rule embed_panel_restricted_genes_all:
    input:
        out_files_panel