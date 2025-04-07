"""
Utility functions during run time.
"""

import os
import re
from typing import Any

COUNT_CORRECTION_MATHOD_WITH_ANNOTATION_PAT = r"(?!(ovrlpy|resolvi_unsupervised)).*"


def uniquify_elements_in_list(
    *x: list[list[Any]],
) -> list[Any]:
    """Remove duplicates in a list of lists.

    Args:
        x (list[list[Any]]): A list of lists.

    Returns:
        list[Any]: A list without duplicates.
    """
    vals: list[Any] = []
    for i in x:
        vals.extend(i)

    return list(dict.fromkeys(vals))


def cross_values_by_key(
    dict_1: dict[str, list[str]],
    dict_2: dict[str, list[str]],
) -> list[tuple[str, str]]:
    """For each common key, cross values in two dictionaries.

    Args:
        dict_1 (dict[str, list[str]]): Left dictionary of values organized by key.
        dict_2 (dict[str, list[str]]): Right dictionary of values organized by key.

    Returns:
        list[list[str]]: Results of crossed dict_1 and dict_2.
    """
    common_keys: set[str] = set(dict_1.keys()) & set(dict_2.keys())

    ret: list[tuple[str, str]] = []

    for k in common_keys:
        if dict_1[k] is None or dict_2[k] is None:
            continue

        assert isinstance(dict_1[k], list) and isinstance(dict_2[k], list)
        if len(dict_1[k]) == 0 or len(dict_2[k]) == 0:
            continue

        ret.extend(
            [(i, j) for i in dict_1[k] for j in dict_2[k]],
        )

    return ret


def get_valid_segmentation_methods_per_count_correction_method(
    candidate_segmentation_methods: list[str],
    count_correction_method: str,
) -> list[str]:
    if re.match(
        r"^ovrlpy$",
        count_correction_method,
        flags=re.IGNORECASE,
    ):
        return [
            i
            for i in candidate_segmentation_methods
            if re.match(
                r"(^10x_\w*?_?0um$)|(^proseg_expected$)",
                i,
                flags=re.IGNORECASE,
            )
            is not None
        ]

    return candidate_segmentation_methods


def get_size(
    path: str,
) -> int:
    """Get the size of a file or directory.

    Args:
        path (str): Path to the file or directory.

    Returns:
        int: Size of the file or directory.
    """
    if not os.path.exists(path):
        raise RuntimeError(f"Error! Input path does not exist: {path}")

    if os.path.isfile(path):
        return os.path.getsize(path)

    total_size: int = 0
    for dirpath, _, filenames in os.walk(path):
        for f in filenames:
            fp: str = os.path.join(dirpath, f)

            if not os.path.islink(fp):
                total_size += os.path.getsize(fp)

    return total_size
