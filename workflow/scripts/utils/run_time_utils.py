"""
Utility functions during run time.
"""


def cross_samples_annotation_by_disease(
    samples: dict[str, list[str]], annotation: dict[str, list[str]]
) -> list[list[str]]:
    """For each disease, apply every annotation method to every sample.

    Args:
        samples (dict[str, list[str]]): A dictionary of samples organized by diseases.
        annotation (dict[str, list[str]]): A dictionary of annotation methods organized by diseases.

    Returns:
        list[list[str]]: Results of crossed samples and annotation.
    """
    assert len(samples) == len(annotation)

    ret: list[list[str]] = []

    for k, v in samples.items():
        assert k in annotation

        ret.extend([[i, j] for i in v for j in annotation[k]])

    return ret
