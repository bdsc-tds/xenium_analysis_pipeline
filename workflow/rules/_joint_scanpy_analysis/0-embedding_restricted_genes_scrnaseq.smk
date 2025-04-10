import pandas as pd

# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
results_dir = Path(config['results_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
scrnaseq_processed_data_dir = Path(config['scrnaseq_processed_data_dir'])
seurat_to_h5_dir = results_dir / 'seurat_to_h5'

# params from pipeline config
min_counts = 10
min_features = 5
max_counts = float("inf")
max_features = float("inf")
min_cells = 5

# Params
layer = 'RNA_counts'
# genes = pd.read_csv('/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/Xenium_hLung_v1_metadata.csv')['Gene'].tolist()
genes = pd.read_csv('/work/PRTNR/CHUV/DIR/rgottar1/spatial/env/xenium_paper/data/markers/Xenium_NSCLC_5k_lung_chromium_common_genes.csv')['gene'].tolist()
n_comps = config['umap_n_comps']
n_neighbors = config['umap_n_neighbors']
min_dist = config['umap_min_dist']
metric = config['umap_metric']


out_files = []

for reference in (references := scrnaseq_processed_data_dir.iterdir()):
    reference_name = reference.stem
    reference_dir = seurat_to_h5_dir / reference_name
    reference_is_done = reference_dir / '.done'

    out_file = results_dir / f'embed_panel_restricted_genes_scrnaseq/{reference_name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet' 
    out_files.append(out_file)

    rule:
        name: f'embed_panel_restricted_genes_scrnaseq/{reference_name}'
        input:
            reference_is_done=reference_is_done
        output:
            out_file=out_file,
        params:
            reference=reference_dir,
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
            mem='100GB',
            runtime='8h',
            # slurm_partition = "gpu",
            # slurm_extra = '--gres=gpu:1',
        conda:
            "spatial"
        shell:
            """
            mkdir -p "$(dirname {output.out_file})"

            python -u workflow/scripts/scRNAseq/embed_panel_scrnaseq.py \
                --reference {params.reference} \
                --out_file {output.out_file} \
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


rule embed_panel_restricted_genes_scrnaseq_all:
    input:
        out_files