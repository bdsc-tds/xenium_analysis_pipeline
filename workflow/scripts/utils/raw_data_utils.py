"""
Utility functions for processing raw data.
"""

from typing import Any

import os
import re
import config_utils
import config_constants as cc


def normalise_path(
    dir_name: str,
    *,
    candidate_paths: tuple[str, ...] = (".", "outs"),
    pat_anchor_file: str = r".*analysis_summary.*\.html",
    pat_flags: re.RegexFlag = re.NOFLAG,
    return_dir: bool = True,
    check_exist: bool = True,
) -> str:
    """Determine the path to a directory based on the specified `anchor_file`.

    Args:
        dir_name (str): Input path to teh raw data directory.
        candidate_paths (tuple[str, ...], optional): Paths to search for `anchor_file`, relative to `dir_name`. Defaults to (".", "outs").
        pat_anchor_file (str, optional): The pattern to a file residing in the valid directory of the raw data. Defaults to ".*analysis_summary.*\\.html".
        pat_flags (re.RegexFlag, optional): Regex flags for the provided pattern. Defaults to re.NOFLAG.
        return_dir (bool, optional): Whether to return the parent directory of the anchor file (`True`) or the path to the file itself (`False`). Defaults to True.

    Raises:
        FileNotFoundError: No files found for the given pattern in the given directory.

    Returns:
        str: The parent directory of the anchor file or the path to the file itself.
    """
    _candidate_dirs: list[str] = [
        os.path.normpath(os.path.join(dir_name, i)) for i in candidate_paths
    ]

    if not check_exist:
        assert len(_candidate_dirs) == 1

        return (
            _candidate_dirs[0]
            if return_dir
            else os.path.normpath(os.path.join(_candidate_dirs[0], pat_anchor_file))
        )

    for i in _candidate_dirs:
        for _, _, filenames in os.walk(i):
            for j in filenames:
                if re.match(pat_anchor_file, j, flags=pat_flags):
                    return i if return_dir else os.path.normpath(os.path.join(i, j))
            break

    raise FileNotFoundError(
        f"Error! The provided file pattern {pat_anchor_file} does not match any files in the following directory: {','.join(_candidate_dirs)}"
    )


def get_gene_panel_file(sample: str, config: dict[str, Any]) -> str | None:
    """Get the path to the provided gene panel file.

    Args:
        sample (str): The sample from which the gene panel file to be extracted.
        config (dict[str, Any]): Configuration dictionary.

    Returns:
        str | None: The path to the gene panel file if provided; otherwise `None` is returned.
    """
    ret: str = config_utils.get_dict_value(
        config,
        "experiments",
        cc.EXPERIMENTS_GENE_PANEL_FILES_NAME,
        config_utils.extract_layers_from_experiments(sample, [0, 1])[0],
    )

    if ret is None:
        return None

    return _get_gene_panel_file(ret)


def _get_gene_panel_file(file_name: str) -> str:
    if os.path.isfile(file_name):
        return file_name

    return os.path.sep.join(
        [
            "/opt/xeniumranger-xenium/lib/json/definitions/panel_designer/panels",
            file_name,
        ]
    )
