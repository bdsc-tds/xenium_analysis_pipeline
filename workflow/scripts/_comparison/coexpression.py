import numpy as np
import pandas as pd
import scipy


def counts_to_positivity(X, cutoff):
    """
    Convert counts to a binary positivity matrix based on a cutoff value.

    Parameters:
    X (scipy.sparse matrix): Input matrix with counts.
    cutoff (float): Cutoff value to determine positivity.

    Returns:
    tuple: A tuple containing the binary positivity matrix and positivity rate.
    """
    positivity = (X >= cutoff).astype(float)
    positivity_rate = positivity.mean(axis=0).A1
    # positivity_rate = np.clip(positivity.mean(axis=0).A1, min_pos_rate, None)
    return positivity, positivity_rate


def conditional_coexpression(positivity):
    """
    Compute conditional co-expression from a positivity matrix.

    Parameters:
    positivity (numpy.ndarray): Binary positivity matrix.

    Returns:
    numpy.ndarray: Conditional co-expression matrix.
    """
    coex = positivity.T @ positivity
    cond_coex = coex / positivity.sum(axis=0)
    return cond_coex.toarray()


def jaccard_coexpression(X):
    """
    Compute the Jaccard index between the rows of X.
    Parameters:
    X (scipy.sparse matrix): Input binary matrix.

    Returns:
    numpy.ndarray: Jaccard index matrix.
    """

    X = X.astype(bool).astype(int)
    intrsct = X.dot(X.T)
    row_sums = intrsct.diagonal()
    unions = row_sums[:, None] + row_sums - intrsct
    jaccard_index = intrsct / unions

    # if min_jaccard > 0.0:
    #     min_jaccard = np.nan_to_num(min_jaccard)
    #     min_jaccard = np.clip(jaccard_index, min_jaccard, None)
    return jaccard_index


def pearson_coexpression(X):
    """
    Compute the Pearson correlation coefficient matrix.

    Parameters:
    X (scipy.sparse matrix): Input matrix.

    Returns:
    numpy.ndarray: Pearson correlation coefficient matrix.
    """
    return np.corrcoef(X.toarray().T)


def spearman_coexpression(X):
    """
    Compute the Spearman rank correlation coefficient matrix.

    Parameters:
    X (scipy.sparse matrix): Input matrix.

    Returns:
    numpy.ndarray: Spearman correlation coefficient matrix.
    """
    return scipy.stats.spearmanr(X.toarray()).statistic

def thin_counts(X, target_count, gen=None):
    """
    Downsample counts to a target count per row.

    Parameters:
    X (scipy.sparse matrix): Input count matrix.
    target_count (int): Target count for each row.
    gen (numpy.random.Generator, optional): Random generator instance.

    Returns:
    scipy.sparse.csr_matrix: Downsampled count matrix.
    """
    if gen is None:
        gen = np.random.default_rng()

    n_counts = X.sum(axis=1)
    probabilities = (X / n_counts).toarray()
    X_thin = np.random.default_rng().multinomial(target_count, probabilities)
    return scipy.sparse.csr_matrix(X_thin)


def sparsify(X):
    """
    Convert a dense matrix to a sparse matrix if not already sparse.

    Parameters:
    X (numpy.ndarray or scipy.sparse matrix): Input matrix.

    Returns:
    scipy.sparse.csr_matrix: Sparse matrix.
    """
    if not scipy.sparse.issparse(X):
        return scipy.sparse.csr_matrix(X.astype(float))
    return X.astype(float)


def coexpression(
    adata,
    positivity_cutoff=1,
    min_samples=0,
    target_count=50,
    method="conditional",
    seed=0,
    # min_positivity_rate: float = 0.01,
    # min_cond_coex: float =0.0,
):
    """
    Calculate co-expression matrix using different methods.

    Parameters:
    adata (anndata.AnnData): AnnData object containing gene expression data.
    positivity_cutoff (float): Cutoff for determining positivity.
    min_samples (int): Minimum number of samples required.
    target_count (int): Target count for downsampling.
    method (str): Method for co-expression calculation ("conditional", "jaccard", "pearson", "spearman").
    seed (int): Seed for random number generation.

    Returns:
    tuple: A tuple containing co-expression matrix, downsampled matrix, positivity matrix, positivity rate, and mask.
    """
    gen = np.random.default_rng(seed)

    X = sparsify(adata.X)

    # Apply mask based on target count threshold
    n_counts = X.sum(axis=1).A1

    if target_count is not None:
        mask = n_counts >= target_count
        print(
            sum(mask),
            "/",
            X.shape[0],
            "(",
            round(100 * sum(mask) / X.shape[0], 2),
            "% ) cells reaching the target count",
        )

        # Modify mask based on min_samples threshold
        if sum(mask) < min_samples:
            print(
                f"Less than {min_samples=} reach the target count. Setting to {min_samples}"
            )
            mask = np.argsort(n_counts)[::-1][:min_samples]

        # Downsample counts
        X_downsample = thin_counts(X[mask], target_count, gen=gen)
    else:
        mask = np.ones(X.shape[0], dtype=bool)
        X_downsample = X

    # Convert counts to binary positivity matrix based on threshold
    pos, pos_rate = counts_to_positivity(
        X_downsample,
        positivity_cutoff,
    )  # min_positivity_rate)

    # Calculate conditional co-expression
    if method == "conditional":
        CC = conditional_coexpression(pos)
    elif method == "jaccard":
        CC = jaccard_coexpression(pos.T).toarray()
    elif method == "pearson":
        CC = pearson_coexpression(X_downsample)
    elif method == "spearman":
        CC = spearman_coexpression(X_downsample)

    CC = pd.DataFrame(CC, index=adata.var_names, columns=adata.var_names)
    pos_rate = pd.Series(pos_rate, index=adata.var_names)

    return CC, X_downsample, pos, pos_rate, mask


def censored_ratio(
    CC_ref_seg,
    CC_other_seg,
    pos_rate_ref_seg=None,
    pos_rate_other_seg=None,
    min_positivity_rate=0.0,
    log2=True,
):
    """
    Compute the ratio of co-expression matrices with optional censoring and log transformation.

    Parameters:
    CC_ref_seg (pd.DataFrame): Reference co-expression matrix.
    CC_other_seg (pd.DataFrame): Other co-expression matrix for comparison.
    pos_rate_ref_seg (pd.Series, optional): Positivity rate for reference segmentation.
    pos_rate_other_seg (pd.Series, optional): Positivity rate for other segmentation.
    min_positivity_rate (float): Minimum positivity rate for filtering.
    log2 (bool): Whether to apply log2 transformation.

    Returns:
    pd.DataFrame: Censored and transformed ratio matrix.
    """
    CCdiff = CC_other_seg / CC_ref_seg
    if log2:
        CCdiff = np.log2(CCdiff)

    # exclude stuff that's barely expressed in one or the other
    if min_positivity_rate > 0.0:
        mask = (pos_rate_ref_seg < min_positivity_rate) | (
            pos_rate_other_seg < min_positivity_rate
        )
        # CCdiff[np.ix_(mask, mask)] = 0
        CCdiff.loc[mask, mask] = 0.0

    return CCdiff


def compare_segmentations(
    CC_ref_seg,
    CC_other_seg,
    pos_rate_ref_seg,
    pos_rate_other_seg,
    min_positivity_rate=0.01,
    cc_cutoff=2.0,
    method=None,
    log2=True,
):
    """
    Compare two co-expression segmentations and identify spurious gene pairs.
    Parameters:
    CC_ref_seg (pd.DataFrame): Reference co-expression matrix.
    CC_other_seg (pd.DataFrame): Other co-expression matrix for comparison.
    pos_rate_ref_seg (pd.Series): Positivity rate for reference segmentation.
    pos_rate_other_seg (pd.Series): Positivity rate for other segmentation.
    min_positivity_rate (float): Minimum positivity rate for filtering.
    cc_cutoff (float): Cutoff for spurious gene pair identification.
    method (str, optional): Method for co-expression calculation.
    log2 (bool): Whether to apply log2 transformation.

    Returns:
    tuple: A tuple containing the difference matrix and spurious gene pairs.
    """

    CCdiff = censored_ratio(
        CC_ref_seg,
        CC_other_seg,
        pos_rate_ref_seg=pos_rate_ref_seg,
        pos_rate_other_seg=pos_rate_other_seg,
        min_positivity_rate=min_positivity_rate,
        log2=log2,
    )

    if log2:
        cc_cutoff = np.log2(cc_cutoff)

    if method == "conditional":
        spurious_gene_pairs = np.where(CCdiff >= cc_cutoff)
    else:
        CCdiff_triu = np.triu(CCdiff, 1)
        spurious_gene_pairs = np.where(CCdiff_triu >= cc_cutoff)

    spurious_gene_pairs = np.array(
        (
            CCdiff.index[spurious_gene_pairs[0]],
            CCdiff.columns[spurious_gene_pairs[1]],
        )
    ).T
    return CCdiff, spurious_gene_pairs
