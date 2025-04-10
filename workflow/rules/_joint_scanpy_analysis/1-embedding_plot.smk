# cfg paths
xenium_dir = Path(config['xenium_processed_data_dir'])
std_seurat_analysis_dir = Path(config['xenium_std_seurat_analysis_dir'])
results_dir = Path(config['results_dir'])
figures_dir = Path(config['figures_dir'])
palette_dir = Path(config['xenium_metadata_dir'])
cell_type_annotation_dir = Path(config['xenium_cell_type_annotation_dir'])

# Params
n_comps = config['umap_n_comps']
n_neighbors = config['umap_n_neighbors']
min_dist = config['umap_min_dist']
metric = config['umap_metric']

s=0.5
alpha=0.5
dpi = 300
points_only = True

cell_type_palette = palette_dir / 'col_palette_cell_types_combo.csv'
panel_palette = palette_dir / 'col_palette_panel.csv'
sample_palette = palette_dir / 'col_palette_sample.csv'

normalisations = ['lognorm','sctransform']
layers = ['data','scale_data']
references = ['matched_reference_combo','external_reference']
methods = ['rctd_class_aware']
colors = ['sample','Level2.1','Level2']#['Level1','Level2','Level3','Level4','panel','sample',] # condition and sample as color to plot added here in addition to levels
extension = 'png'
plot_types = ['umap','facet_umap']

out_files_panel = []

for plot_type in plot_types:
    for segmentation in (segmentations := std_seurat_analysis_dir.iterdir()):
        if segmentation.stem == 'proseg_mode':
            continue
        for condition in (conditions := segmentation.iterdir()): 
            for panel in (panels := condition.iterdir()):
                for normalisation in normalisations:
                    for layer in layers: 
                        for reference in references:
                            for method in methods:
                                for color in colors:
                                    if color == 'Level2.1' and reference == 'external_reference':
                                        continue
                                    if color == 'Level2' and reference == 'matched_reference_combo':
                                        continue

                                    # input embedding file (doesn't depend on ref,method or color loops but more readable to have here)
                                    k = (segmentation.stem,condition.stem,panel.stem,normalisation)
                                    name = '/'.join(k)
                                    embed_file = results_dir / f'embed_panel/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet'

                                    # no need to plot panel for panel color UMAPs
                                    if color == 'panel':
                                        continue
                                    
                                    # no need to plot sample coloring for every param combination or for facet plots
                                    if color == 'sample':
                                        if (reference != references[0] or method != methods[0]):
                                            continue
                                        if plot_type == 'facet_umap':
                                            continue


                                    out_file = figures_dir / f"embed_panel/{name}/{plot_type}_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{reference}_{method}_{color}.{extension}"
                                    out_files_panel.append(out_file)

                                    rule:
                                        name: f'embed_panel_plot/{name}/{plot_type}_{layer}_{reference}_{method}_{color}'
                                        input:
                                            panel=panel,
                                            embed_file=embed_file,
                                        output:
                                            out_file=out_file,
                                        params:
                                            cell_type_annotation_dir=cell_type_annotation_dir,
                                            normalisation=normalisation,
                                            reference=reference,
                                            method=method,
                                            color=color,
                                            cell_type_palette=cell_type_palette,
                                            panel_palette=panel_palette,
                                            sample_palette=sample_palette,
                                            s=s,
                                            alpha=alpha,
                                            dpi=dpi,
                                            points_only='--points_only' if points_only else '',
                                            facet='--facet' if plot_type == 'facet_umap' else ''
                                        threads: 1
                                        resources:
                                            mem='30GB',
                                            runtime='10m',
                                        conda:
                                            "spatial"
                                        shell:
                                            """
                                            mkdir -p "$(dirname {output.out_file})"

                                            python workflow/scripts/xenium/embed_panel_plot.py \
                                            --panel {input.panel} \
                                            --embed_file {input.embed_file} \
                                            --cell_type_annotation_dir {params.cell_type_annotation_dir} \
                                            --normalisation {params.normalisation} \
                                            --reference {params.reference} \
                                            --method {params.method} \
                                            --color {params.color} \
                                            --out_file {output.out_file} \
                                            --cell_type_palette {params.cell_type_palette} \
                                            --panel_palette {params.panel_palette} \
                                            --sample_palette {params.sample_palette} \
                                            --s {params.s} \
                                            --alpha {params.alpha} \
                                            --dpi {params.dpi} \
                                            {params.points_only} \
                                            {params.facet} \
                                            
                                            echo "DONE"
                                            """






# out_files_condition = []
# for segmentation in (segmentations := xenium_dir.iterdir()):
#     if segmentation.stem == 'proseg_v1':
#         continue
#     for condition in (conditions := segmentation.iterdir()): 
#         for panel in (panels := condition.iterdir()):
#             for normalisation in normalisations:
#                 for layer in layers: 

#                     k = (segmentation.stem,condition.stem,normalisation)
#                     name = '/'.join(k)
#                     embed_file = results_dir / f'embed_condition/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}.parquet'

#                     for reference in references:
#                         for method in methods:
#                             for color in colors:
                                
#                                 # no need to plot sample coloring for every param combination
#                                 if color == 'sample' and reference != references[0] and method != methods[0]:
#                                     continue

#                                 out_file = figures_dir / f"embed_condition/{name}/umap_{layer}_{n_comps=}_{n_neighbors=}_{min_dist=}_{metric}_{reference}_{method}_{color}.{extension}"
#                                 out_files_condition.append(out_file)

#                                 rule:
#                                     name: f'embed_condition_plot/{name}/umap_{layer}_{reference}_{method}_{color}'
#                                     input:
#                                         condition=condition,
#                                         embed_file=embed_file,
#                                     output:
#                                         out_file=out_file,
#                                     params:
#                                         cell_type_annotation_dir=cell_type_annotation_dir,
#                                         reference=reference,
#                                         method=method,
#                                         color=color,
#                                         cell_type_palette=cell_type_palette,
#                                         panel_palette=panel_palette,
#                                         sample_palette=sample_palette,
#                                         s=s,
#                                         alpha=alpha,
#                                     threads: 1
#                                     resources:
#                                         mem='30GB',
#                                         runtime='10m',
#                                     conda:
#                                         "spatial"
#                                     shell:
#                                         """
#                                         mkdir -p "$(dirname {output.out_file})"

#                                         python workflow/scripts/xenium/embed_condition_plot.py \
#                                         --condition {input.condition} \
#                                         --embed_file {input.embed_file} \
#                                         --cell_type_annotation_dir {params.cell_type_annotation_dir} \
#                                         --reference {params.reference} \
#                                         --method {params.method} \
#                                         --color {params.color} \
#                                         --out_file {output.out_file} \
#                                         --cell_type_palette {params.cell_type_palette} \
#                                         --panel_palette {params.panel_palette} \
#                                         --sample_palette {params.sample_palette} \
#                                         --s {params.s} \
#                                         --alpha {params.alpha} \

#                                         echo "DONE"
#                                         """





rule embed_panel_plot_all:
    input:
        out_files_panel

# rule embed_condition_plot_all:
#     input:
#         out_files_condition