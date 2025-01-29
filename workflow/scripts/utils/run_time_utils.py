"""
Utility functions during run time.
"""

from typing import Any


def uniquify_elements_in_list(*x: list[list[Any]]) -> list[Any]:
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
    vals_1: dict[str, list[str]], vals_2: dict[str, list[str]]
) -> list[tuple[str, str]]:
    """For each key, cross their values.

    Args:
        vals_1 (dict[str, list[str]]): Left dictionary of values organized by key.
        vals_2 (dict[str, list[str]]): Right dictionary of values organized by key.

    Returns:
        list[list[str]]: Results of crossed vals_1 and vals_2.
    """
    assert len(vals_1) == len(vals_2)

    ret: list[tuple[str, str]] = []

    for k, v in vals_1.items():
        assert k in vals_2

        ret.extend([(i, j) for i in v for j in vals_2[k]])

    return ret
