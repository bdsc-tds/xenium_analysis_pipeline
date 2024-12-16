"""
Utility functions for processing the configuration dictionary to prepare for the execution of the workflow. The configuration dictionary parsed automatically by Snakemake will be modified in place.
"""

from typing import Any, Callable

import re
import os
import yaml

import config_constants as cc


def extract_layers_from_experiments(
    experiments: list[str] | str,
    layers: list[int] | int,
    sep_in: str = os.path.sep,
    sep_out: str | None = os.path.sep,
) -> list[Any]:
    """Extract layers from experiments.

    Args:
        experiments (list[str] | str): Input experiments.
        layers (list[int] | int): Layers to be extracted.
        sep_in (str, optional): Separator used in the input experiments. Defaults to os.path.sep.
        sep_out (str | None, optional): Separator to be used for the extracted layers. Defaults to os.path.sep.

    Returns:
        list[Any]: A list of extracted layers.
    """
    ret: list[Any] = []

    experiments = _convert2list(experiments, match_length=False)
    layers = sorted(_convert2list(layers, match_length=False))

    for i in experiments:
        j = i.split(sep_in)
        assert max(layers) < len(j)

        if sep_out is not None:
            ret.append(sep_out.join([j[k] for k in layers]))
        else:
            ret.append(j)

    return ret


def set_dict_value(
    x: dict[str | int | float | tuple, Any],
    *keys: str | int | float | tuple,
    value: Any,
    force: bool = True,
):
    """Assign the value in a dictionary in place located by keys.

    Args:
        x (dict[str  |  int  |  float  |  tuple, Any]): A dictionary.
        *keys (str | int | float | tuple): Keys in the given dictionary to locate the value.
        value (Any): The value to be assigned.
        force (bool, optional): Whether to create an empty dictionary if certain keys are not found. Defaults to True.

    Raises:
        KeyError: If some keys are not found in the given dictionary and an empty dictionary will not be created.
    """
    for k in range(len(keys) - 1):
        if keys[k] not in x:
            if force:
                x[keys[k]] = {}
            else:
                raise KeyError(f"Error! {keys[k]} is not in {list(x.keys())}")

        x = x[keys[k]]

    x[keys[-1]] = value


def get_dict_value(
    x: dict[Any, Any] | None,
    *keys: str | int | float | tuple,
    replace_none: Any = None,
    inexist_key_ok: bool = False,
) -> Any:
    """Get a value from a given dictionary using given keys.

    Args:
        x (dict[Any, Any]): A dictionary from which a value is retrieved.
        replace_none (Any, optional): Replacement value if the retrieved value is `None`. Defaults to None.

    Raises:
        KeyError: When a given key is not in the dictionary.

    Returns:
        Any: Retrieved value.
    """
    if x is not None:
        for k in keys:
            if k not in x:
                if inexist_key_ok:
                    x = None
                    break
                raise KeyError(f"Error! Key {k} is not in {list(x.keys())}")

            x = x[k]

    if x is None:
        return replace_none
    return x


def _convert2list(
    x: Any, length: int = 1, *, match_length: bool = True, convert2list: bool = True
) -> list | tuple | set:
    assert length > 0

    if not isinstance(x, (list, set, tuple)) or not match_length and length > 1:
        return [x] * length

    if match_length and len(x) == 1:
        return (x if isinstance(x, list) else list(x)) * length

    if (not match_length and length == 1) or (match_length and len(x) == length):
        return x if (isinstance(x, list) or not convert2list) else list(x)

    raise RuntimeError(f"Error! Cannot convert {x} to a list with length {length}.")


def _merge_flattened_lists(
    left: list,
    right: list,
    *,
    pad_layers: dict[int, int | float | str | None] | None = None,
) -> list:
    if len(left) < len(right):
        less = left
        more = right
    else:
        less = right
        more = left

    if len(less) == 0:
        return more

    ret: list = []

    if pad_layers is None:
        assert len(left) == len(right)
    else:
        if len(left) == len(right):
            pad_layers = None

    idx_less: int = 0
    for idx, e in enumerate(more):
        assert isinstance(e, list)

        if pad_layers is None or idx not in pad_layers:
            assert isinstance(less[idx_less], list)

            ret.append([*e, *less[idx_less]])
            idx_less += 1
        else:
            ret.append(
                [
                    *e,
                    *_convert2list(pad_layers[idx], len(less[0])),
                ]
            )

    return ret


def _merge_dicts(
    x: list[dict | None] | set[dict | None] | tuple[dict | None], exist_ok: bool = False
) -> dict | None:
    def __merge_to(keep: Any, add: Any, *, exist_ok: bool = False) -> dict | list:
        assert keep is not None and add is not None

        if add is None or isinstance(add, (dict, list)) and len(add) == 0:
            return keep

        if isinstance(keep, dict) + isinstance(add, dict) == 1:
            raise RuntimeError(
                f"Error! Cannot merge as at least one of the following is not a dictionary: {keep}, {add}"
            )

        if isinstance(keep, dict) + isinstance(add, dict) == 0:
            return [
                *_convert2list(keep, match_length=False),
                *_convert2list(add, match_length=False),
            ]

        for k, v in add.items():
            if k is None:
                continue

            if k in keep:
                if exist_ok:
                    keep[k] = __merge_to(keep[k], v, exist_ok=exist_ok)
                else:
                    raise RuntimeError(f"Error! Duplicate keys {k} found.")
            else:
                keep[k] = v

        return keep

    if len(x) == 0:
        return None

    ret: dict = {}

    for i in x:
        if i is None:
            continue

        ret = __merge_to(ret, i, exist_ok=exist_ok)

    return ret if len(ret) > 0 else None


def _flatten_struct(
    x: Any,
    *,
    layer: int = 0,
    key_layer2pat_include: dict[tuple[int, ...] | str, str | Callable] | None = None,
    key_layer2pat_exclude: dict[tuple[int, ...] | str, str | Callable] | None = None,
    pad_layers: dict[int, int | float | str | None] | None = None,
) -> tuple[list[Any], ...]:
    if key_layer2pat_include is None:
        key_layer2pat_include = {"_": r"^(?!_+).+"}

    def __is_key4layer(
        _key: str,
        _layer: int,
        key_layer2pat: dict[tuple[int, ...] | str, str | Callable],
    ) -> bool:
        assert key_layer2pat is not None
        for k, v in key_layer2pat.items():
            if k == "_" or _layer in k:
                if isinstance(v, str):
                    return re.match(v, _key) is not None
                return v(_key)
        return False

    ret: list = []

    if isinstance(x, dict):
        for key, value in x.items():
            if key_layer2pat_exclude is not None and __is_key4layer(
                key, layer, key_layer2pat_exclude
            ):
                continue

            if __is_key4layer(key, layer, key_layer2pat_include):
                flattened = _flatten_struct(
                    value,
                    layer=layer + 1,
                    key_layer2pat_include=key_layer2pat_include,
                    key_layer2pat_exclude=key_layer2pat_exclude,
                    pad_layers=pad_layers,
                )
                flattened = _convert2list(key, len(flattened[0])), *list(flattened)

                ret = _merge_flattened_lists(
                    ret, list(flattened), pad_layers=pad_layers
                )

        return tuple(ret)
    else:
        _ret = _convert2list(x, len(x) if isinstance(x, (list, set, tuple)) else 1)

        for i in _ret:
            if (
                i is not None
                and key_layer2pat_exclude is not None
                and __is_key4layer(i, layer, key_layer2pat_exclude)
            ):
                continue

            if i is None or __is_key4layer(i, layer, key_layer2pat_include):
                ret.append(i)

        return (ret,)


def _collect_experiments(
    flattened: tuple[list[str | int | float], ...],
    layer: int,
    *,
    drop_layers: tuple[int, ...] | None = None,
    key_func: Callable = os.path.sep.join,
    val_func: Callable | None = os.path.sep.join,
    key_val_func: Callable | None = lambda ks, vs: os.path.sep.join(ks + vs),
    vals_func: Callable | None = None,
    simplify: Callable | None = None,
) -> dict[str, Any]:
    assert 0 <= layer < len(flattened)

    for i in range(1, len(flattened)):
        assert len(flattened[0]) == len(flattened[i])

    ret: dict[str, Any] = {}

    for idx, _ in enumerate(flattened[0]):
        k: str = key_func(
            [
                flattened[i][idx]
                for i in range(layer + 1)
                if drop_layers is None or i not in drop_layers
            ]
        )

        if k == "":
            raise KeyError("Error! Invalid empty key generated.")

        v: tuple[str | int | float, ...] | str | int | float = ()
        for i in range(layer + 1, len(flattened)):
            if drop_layers is None or i not in drop_layers:
                v += (flattened[i][idx],)

        if val_func is not None:
            v = val_func(v)

        if key_val_func is not None:
            v = key_val_func(
                _convert2list(k, match_length=False),
                _convert2list(v, match_length=False, convert2list=False),
            )

        if isinstance(v, (list, tuple)) and len(v) == 1:
            v = v[0]

        if k not in ret:
            ret[k] = []
        ret[k].append(v)

    if vals_func is not None:
        for key, value in ret.items():
            ret[key] = vals_func(value)

    if simplify is not None:
        return simplify(ret)

    return ret


def _process_experiments(file_path: str, root_path: str) -> tuple[Any, ...]:
    if not os.path.isabs(file_path):
        file_path = os.path.join(root_path, file_path)

    if not os.path.isfile(file_path):
        raise FileNotFoundError(
            f"Error! File containing experiment information does not exist: {file_path}"
        )

    with open(file_path, "r", encoding="utf-8") as fh:
        data = yaml.safe_load(fh)

    flattened = _flatten_struct(
        data,
        key_layer2pat_exclude={
            "_": "|".join(
                [
                    cc.EXPERIMENTS_CONFIG_PATH_NAME,
                    cc.EXPERIMENTS_BASE_PATH_NAME,
                    cc.EXPERIMENTS_COLLECTIONS_NAME,
                    cc.EXPERIMENTS_GENE_PANEL_FILES_NAME,
                    cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
                    cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
                ]
            )
        },
    )

    diseases_level: dict[str, Any] = _collect_experiments(flattened, 0)
    gene_panels_level: dict[str, Any] = _collect_experiments(flattened, 1)
    donors_level: dict[str, Any] = _collect_experiments(flattened, 2)
    samples_level: dict[str, Any] = _collect_experiments(flattened, 3)

    wildcards: dict[str, Any] = {}
    wildcards[cc.WILDCARDS_DISEASES_NAME] = list(diseases_level.keys())
    wildcards[cc.WILDCARDS_GENE_PANELS_NAME] = list(gene_panels_level.keys())
    wildcards[cc.WILDCARDS_DONORS_NAME] = list(donors_level.keys())
    wildcards[cc.WILDCARDS_SAMPLES_NAME] = list(samples_level.keys())

    collections: dict[str, Any] = {}
    collections[cc.EXPERIMENTS_COLLECTIONS_DISEASES_NAME] = diseases_level
    collections[cc.EXPERIMENTS_COLLECTIONS_GENE_PANELS_NAME] = gene_panels_level
    collections[cc.EXPERIMENTS_COLLECTIONS_DONORS_NAME] = donors_level

    gene_panel_files: dict[str, Any] = _collect_experiments(
        _flatten_struct(
            data,
            key_layer2pat_include={
                (0, 1, 3): r"^(?!_+).+",
                (2,): cc.EXPERIMENTS_GENE_PANEL_FILES_NAME,
            },
            key_layer2pat_exclude={
                "_": "|".join(
                    [
                        cc.EXPERIMENTS_CONFIG_PATH_NAME,
                        cc.EXPERIMENTS_BASE_PATH_NAME,
                        cc.EXPERIMENTS_COLLECTIONS_NAME,
                        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
                        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
                    ]
                )
            },
        ),
        1,
        drop_layers=(2,),
        val_func=None,
        key_val_func=None,
        simplify=lambda x: {k: v[0] if len(v) == 1 else v for k, v in x.items()},
    )

    gene_panel_qc_thresholds = _collect_experiments(
        _flatten_struct(
            data,
            key_layer2pat_include={
                (0, 1, 3): r"^(?!_+).+",
                (2,): cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
                (4,): lambda v: (
                    v == "Inf" if isinstance(v, str) else isinstance(v, (int, float))
                ),
            },
            key_layer2pat_exclude={
                (0, 1, 2, 3): "|".join(
                    [
                        cc.EXPERIMENTS_CONFIG_PATH_NAME,
                        cc.EXPERIMENTS_BASE_PATH_NAME,
                        cc.EXPERIMENTS_COLLECTIONS_NAME,
                        cc.EXPERIMENTS_GENE_PANEL_FILES_NAME,
                        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
                    ]
                ),
            },
            pad_layers={3: None},
        ),
        1,
        drop_layers=(2,),
        val_func=lambda vs: (
            None if vs[0] is None else ({vs[0]: vs[1]} if len(vs) > 1 else vs[0])
        ),
        key_val_func=None,
        vals_func=_merge_dicts,
        simplify=lambda x: {k: v for k, v in x.items() if v is not None and len(v) > 0},
    )

    cell_type_annotation = _collect_experiments(
        _flatten_struct(
            data,
            key_layer2pat_include={
                (0, 2, 3, 4): r"^(?!_+).+",
                (1,): cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
                (5,): lambda v: isinstance(v, (int, str)),
            },
            key_layer2pat_exclude={
                (0, 1, 2, 3, 4): "|".join(
                    [
                        cc.EXPERIMENTS_CONFIG_PATH_NAME,
                        cc.EXPERIMENTS_BASE_PATH_NAME,
                        cc.EXPERIMENTS_COLLECTIONS_NAME,
                        cc.EXPERIMENTS_GENE_PANEL_FILES_NAME,
                        cc.EXPERIMENTS_GENE_PANEL_QC_NAME,
                    ]
                ),
            },
        ),
        0,
        drop_layers=(1,),
        val_func=lambda v: (
            {
                os.path.join(v[0], v[1]): {
                    v[2]: str(v[3]) if isinstance(v[3], int) else v[3]
                }
            }
        ),
        key_val_func=None,
        vals_func=lambda x: _merge_dicts(x, exist_ok=True),
    )

    return (
        wildcards,
        data[cc.EXPERIMENTS_BASE_PATH_NAME],
        collections,
        gene_panel_files,
        gene_panel_qc_thresholds,
        cell_type_annotation,
    )


def _process_segmentation(data: dict[str, Any]) -> tuple[list[str], dict[str, Any]]:
    methods: list[str] = []
    ret: dict[str, Any] = {}

    _methods: list[str] = get_dict_value(data, "methods")

    if len(_methods) == 0:
        raise RuntimeError("Error! At least one segmentation method should be used.")

    for m in _methods:
        if m == "10x":
            exp_vals: list[int] = _convert2list(
                get_dict_value(data, "10x", "expansion-distance"), match_length=False
            )

            if len(exp_vals) == 0:
                raise RuntimeError(
                    "Error! At least one value for 'expansion-distance' in the config file must be specified."
                )

            methods_10x = ["_".join(["10x", "".join([str(i), "um"])]) for i in exp_vals]
            methods.extend(methods_10x)

            tmp = {
                k: _convert2list(v, len(exp_vals))
                for k, v in get_dict_value(data, "10x").items()
            }

            for idx, val in enumerate(methods_10x):
                ret[val] = {}

                for k, v in tmp.items():
                    ret[val][k] = v[idx]
        else:
            methods.append(m)

    return methods, ret


def _process_cell_type_annotation(
    data_from_experiments: dict[str, Any], data_from_config: dict[str, Any]
) -> tuple[dict[str, dict[str, str]], dict[str, dict[str, str]], dict[str, list[str]]]:
    paths: dict[str, dict[str, str]] = {}
    levels: dict[str, dict[str, str]] = {}
    wildcards: dict[str, list[str]] = {}

    for k_1, v_1 in data_from_experiments.items():
        for k_2, v_2 in list(v_1.items()):
            if "path" not in v_2 or v_2["path"] is None or v_2["path"] == "":
                del data_from_experiments[k_1][k_2]
                continue

            if "levels" not in v_2 or v_2["levels"] is None or len(v_2["levels"]) == 0:
                raise RuntimeError(
                    f"Error! Entry 'levels' of {k_1}'s {k_2} should be specified."
                )

            # Collect paths to reference.
            if k_1 not in paths:
                paths[k_1] = {}
            assert k_2 not in paths[k_1]
            paths[k_1][k_2] = v_2["path"]

            # Collect levels for reference.
            if k_1 not in levels:
                levels[k_1] = {}
            assert k_2 not in levels[k_1]
            levels[k_1][k_2] = v_2["levels"]

            # Generate wildcards.
            approach: str = extract_layers_from_experiments(k_2, 0)[0]

            _wildcards = [
                os.path.join(i[0], j, i[1], m)
                for i in [
                    [k_2, _v] for _v in _convert2list(v_2["levels"], match_length=False)
                ]
                for j in data_from_config[approach]["methods"]
                for m in data_from_config[approach]["modes"]
            ]

            if k_1 not in wildcards:
                wildcards[k_1] = _wildcards
            else:
                wildcards[k_1].extend(_wildcards)

        if len(v_1) == 0:
            raise RuntimeError(
                f"Error! At least one reference type among 'matched_reference' and 'external_reference' should be provided for disease {k_1}."
            )

    for k_1, v_1 in data_from_config.items():
        for k_2, v_2 in v_1.items():
            if k_2 in ["methods", "modes"]:
                continue

            updated_method: dict = {}

            for k_3, v_3 in v_2.items():
                tmp = _convert2list(
                    v_3, len(v_1["modes"]), match_length=len(v_1["modes"]) != 1
                )

                for idx, _v in enumerate(v_1["modes"]):
                    if _v not in updated_method:
                        updated_method[_v] = {}

                    assert k_3 not in updated_method[_v]

                    updated_method[_v][k_3] = tmp[idx]

            set_dict_value(data_from_config, k_1, k_2, value=updated_method)

    return paths, levels, wildcards


def process_config(
    data: dict[str | int | float | tuple, Any], *, root_path: str
) -> None:
    """Process configuration dictionary parsed by Snakemake.

    Processed configuration will be saved in the same dictionary. The detailed steps are as follows:

    1. Process `experiments` section.

    The configuration file containing experiment information will be loaded. Those four layers (diseases, gene panels, donors, and samples) will be flattened under different levels, which can be used as wildcards in the workflow.

    Samples corresponding to the wildcards flattened from the first three layers are available under `_collections` of `experiments`.

    2. Process `segmentation` section.

    Different segmentation methods will be used as wildcards. Particularly, multiple expansion distances can be used for the 10x method (Xenium Ranger), and thus each of them is treated as an independend method for segmentation.

    Args:
        data (dict[str, Any]): Configuration dictionary parsed from an yaml file.
        root_path (str): The root path to the yaml file.
    """

    if cc.WILDCARDS_NAME in data:
        return None

    # Process `experiments` section.
    _experiments = _process_experiments(
        get_dict_value(data, "experiments"),
        root_path,
    )

    set_dict_value(
        data,
        cc.WILDCARDS_NAME,
        value=_experiments[0],
    )

    set_dict_value(
        data,
        "experiments",
        value={
            cc.EXPERIMENTS_CONFIG_PATH_NAME: get_dict_value(data, "experiments"),
            cc.EXPERIMENTS_BASE_PATH_NAME: _experiments[1],
            cc.EXPERIMENTS_COLLECTIONS_NAME: _experiments[2],
            cc.EXPERIMENTS_GENE_PANEL_FILES_NAME: _experiments[3],
            cc.EXPERIMENTS_GENE_PANEL_QC_NAME: _experiments[4],
        },
    )

    # Process `segmentation` section.
    _segmentation = _process_segmentation(get_dict_value(data, "segmentation"))

    set_dict_value(
        data, cc.WILDCARDS_NAME, cc.WILDCARDS_SEGMENTATION_NAME, value=_segmentation[0]
    )

    for k, v in _segmentation[1].items():
        set_dict_value(data, "segmentation", k, value=v)

    # Process `cell_type_annotation` section.
    _cell_type_annotation = _process_cell_type_annotation(
        _experiments[5], get_dict_value(data, "cell_type_annotation")
    )

    set_dict_value(
        data,
        "experiments",
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_REFERENCE_FILES_NAME,
        value=_cell_type_annotation[0],
    )

    set_dict_value(
        data,
        "experiments",
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_NAME,
        cc.EXPERIMENTS_CELL_TYPE_ANNOTATION_LEVELS_NAME,
        value=_cell_type_annotation[1],
    )

    set_dict_value(
        data,
        cc.WILDCARDS_NAME,
        cc.WILDCARDS_CELL_TYPE_ANNOTATION_NAME,
        value=_cell_type_annotation[2],
    )
