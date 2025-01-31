rule summarize_metadata:
    input:
        seurat_objects = expand("segmentation/{condition}/{panel}/{donor}/{sample_id}/cell_type_annotation/{ref_based}/{external}/{class_aware}/{annotation_level}/seurat_object.rds", panel=["panel1", "panel2", "panel3", "panel4"]),
        annotation_file = "annotations/{panel}_annotations.csv"
    output:
        summary_csv = "segmentation/{condition}/{panel}/metadata_summary.csv"
    group: "panel_summary"
    run:
        Rscript("scripts/summarize_metadata.R", input.seurat_objects, input.annotation_file, output.summary_csv)
