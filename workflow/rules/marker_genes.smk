#######################################
#              Functions              #
#######################################

def get_marker_genes_config(*keys: str, default: Any) -> Any:
    """
    Get a value from the `marker_genes` section of the configuration.

    The default is returned if any level is missing or left blank.
    """
    value: Any = config

    for key in ("marker_genes", *keys):
        value = get_dict_value(
            value,
            key,
            inexist_key_ok=True,
        )

    return default if value is None else value


def get_rctd_results4marker_genes(wildcards, prefix: str) -> list[str]:
    """
    Get the path to the RCTD results table, which is used to select singlets.

    An empty list is returned for annotation methods other than RCTD, which have no such table, and when cells are not restricted to singlets, in which case the table is not needed.

    `prefix` is the root directory of cell type annotation results, either before or after count correction. In the latter case the layout of the path depends on the count correction method.
    """
    annotation_method: str = extract_layers_from_experiments(wildcards.annotation_id, [2])[0]

    if re.match(
        r"^rctd_.+",
        annotation_method,
        flags=re.IGNORECASE,
    ) is None:
        return []

    if not get_marker_genes_config("rctd_singlets_only", default=True):
        return []

    suffix: str = "output/results_df.parquet"

    _prefix: str = os.path.join(
        prefix,
        wildcards.segmentation_id,
        wildcards.sample_id,
    )

    count_correction_id: str | None = getattr(wildcards, "count_correction_id", None)

    # Before count correction.
    if count_correction_id is None:
        return [
            os.path.join(
                _prefix,
                wildcards.normalisation_id,
                wildcards.annotation_id,
                suffix,
            )
        ]

    if re.fullmatch(
        r"^ovrlpy$",
        count_correction_id,
        flags=re.IGNORECASE,
    ) is not None:
        return [
            os.path.join(
                _prefix,
                "ovrlpy",
                f'signal_integrity_threshold={get_dict_value(
                    config,
                    "count_correction",
                    "ovrlpy",
                    "signal_integrity_threshold",
                )}',
                wildcards.normalisation_id,
                wildcards.annotation_id,
                suffix,
            )
        ]
    elif re.fullmatch(
        r"^resolvi_unsupervised$",
        count_correction_id,
        flags=re.IGNORECASE,
    ) is not None:
        return [
            os.path.join(
                _prefix,
                "resolvi_unsupervised",
                f'mixture_k={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "train",
                    "mixture_k",
                )}',
                f'num_samples={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "predict",
                    "num_samples",
                )}',
                wildcards.normalisation_id,
                wildcards.annotation_id,
                suffix,
            )
        ]
    elif re.fullmatch(
        r"^resolvi_supervised$",
        count_correction_id,
        flags=re.IGNORECASE,
    ) is not None:
        return [
            os.path.join(
                _prefix,
                wildcards.normalisation_id,
                wildcards.annotation_id,
                "resolvi_supervised",
                f'mixture_k={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "train",
                    "mixture_k",
                )}',
                f'num_samples={get_dict_value(
                    config,
                    "count_correction",
                    "resolvi",
                    "predict",
                    "num_samples",
                )}',
                wildcards.normalisation_id,
                wildcards.annotation_id,
                suffix,
            )
        ]

    # Count correction by SPLIT.
    return [
        os.path.join(
            _prefix,
            wildcards.normalisation_id,
            wildcards.annotation_id,
            count_correction_id,
            wildcards.normalisation_id,
            wildcards.annotation_id,
            suffix,
        )
    ]


######################################
#              Subrules              #
######################################

include: "_marker_genes/scran_score_markers.smk"
