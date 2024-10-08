"""
Utility functions for processing the configuration dictionary to prepare for the execution of the workflow. The configuration dictionary parsed automatically by Snakemake will be modified in place.
"""

from typing import Any

import os
import yaml

import config_constants as cc


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
    x: dict[Any, Any], *keys: str | int | float | tuple, replace_none: Any = None
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

    for k in keys:
        if k not in x:
            raise KeyError(f"Error! Key {k} is not in {list(x.keys())}")

        x = x[k]

    if x is None:
        return replace_none
    return x


def _convert2list(x: Any, length: int = 1, *, match_length: bool = True) -> list[Any]:
    assert length > 0

    if not isinstance(x, (list, set, tuple)) or not match_length and length > 1:
        return [x] * length

    if match_length and len(x) == 1:
        return (x if isinstance(x, list) else list(x)) * length

    if (not match_length and length == 1) or (match_length and len(x) == length):
        return x if isinstance(x, list) else list(x)

    raise RuntimeError(f"Error! Cannot convert {x} to a list with length {length}.")


def _flatten_struct(
    x: Any, *keys2exclude: str | int | float | tuple
) -> tuple[list[Any], ...]:
    ret: list = []

    if isinstance(x, dict):
        for key, value in x.items():
            if key in keys2exclude:
                continue

            flattened = _flatten_struct(value, keys2exclude)
            flattened = _convert2list(key, len(flattened[0])), *list(flattened)

            if len(ret) == 0:
                ret = list(flattened)
            else:
                assert len(ret) == len(flattened)

                for idx, e in enumerate(ret):
                    assert isinstance(e, list) and isinstance(flattened[idx], list)

                    e.extend(flattened[idx])

        return tuple(ret)

    return (_convert2list(x, len(x) if isinstance(x, (list, set, tuple)) else 1),)


def _collect_experiments(
    flattened: tuple[list[str], ...], level: int
) -> dict[str, Any]:
    assert 0 <= level < len(flattened)

    for i in range(1, len(flattened)):
        assert len(flattened[0]) == len(flattened[i])

    ret: dict[str, Any] = {}

    for idx, _ in enumerate(flattened[0]):
        k: str = os.path.join(*[flattened[i][idx] for i in range(level + 1)])
        v: str = os.path.join(
            k, *[flattened[i][idx] for i in range(level + 1, len(flattened))]
        )

        if not k in ret:
            ret[k] = []
        ret[k].append(v)

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
        cc.EXPERIMENTS_CONFIG_PATH_NAME,
        cc.EXPERIMENTS_BASE_PATH_NAME,
        cc.EXPERIMENTS_COLLECTIONS_NAME,
    )

    diseases_level: dict[str, Any] = _collect_experiments(flattened, 0)
    gene_panels_lebel: dict[str, Any] = _collect_experiments(flattened, 1)
    donors_level: dict[str, Any] = _collect_experiments(flattened, 2)
    samples_level: dict[str, Any] = _collect_experiments(flattened, 3)

    wildcards: dict[str, Any] = {}
    wildcards[cc.WILDCARDS_DISEASES_NAME] = list(diseases_level.keys())
    wildcards[cc.WILDCARDS_GENE_PANELS_NAME] = list(gene_panels_lebel.keys())
    wildcards[cc.WILDCARDS_DONORS_NAME] = list(donors_level.keys())
    wildcards[cc.WILDCARDS_SAMPLES_NAME] = list(samples_level.keys())

    collections: dict[str, Any] = {}
    collections[cc.EXPERIMENTS_COLLECTIONS_DISEASES_NAME] = diseases_level
    collections[cc.EXPERIMENTS_COLLECTIONS_GENE_PANELS_NAME] = gene_panels_lebel
    collections[cc.EXPERIMENTS_COLLECTIONS_DONORS_NAME] = donors_level

    return wildcards, data[cc.EXPERIMENTS_BASE_PATH_NAME], collections


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
        data (dict[str, Any]): _description_
        root_path (str): _description_
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
        },
    )

    # Process `segmentation` section.
    _segmentation = _process_segmentation(get_dict_value(data, "segmentation"))

    set_dict_value(
        data, cc.WILDCARDS_NAME, cc.WILDCARDS_SEGMENTATION_NAME, value=_segmentation[0]
    )

    for k, v in _segmentation[1].items():
        set_dict_value(data, "segmentation", k, value=v)

    return None

