from collections import defaultdict
import logging
from pathlib import Path
from typing import AnyStr, Dict, Union, Iterable, Sequence, Tuple, Optional

from anndata import AnnData
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from dataforest.utils.analysis.pairwise import pairwise_metric_heat, pairwise_metric
from scipy.stats import pearsonr

from cellforest.api.io import read_10x
from cellforest.utils.scanpy.ambient import est_ambient_rna
from cellforest.utils.scanpy.doublet import dub_finder_disk
from cellforest.utils.scanpy.group import groupby_dict


def pseudo_bulk(ad: AnnData, groupby: Optional[str] = None, layer="X") -> pd.DataFrame:
    ad_d = groupby_dict(ad, groupby)
    bulk_d = dict()
    for k, _ad in ad_d.items():
        expr = _pseudo_bulk_vector(_ad, layer)
        bulk_d[k] = expr
    df = pd.DataFrame(bulk_d).T
    df.columns = ad.var_names
    return df


def _pseudo_bulk_vector(ad: AnnData, layer="X") -> np.ndarray:
    X = ad.X if layer == "X" else ad.layers[layer]
    expr = X.sum(axis=0)
    expr = np.array(expr / expr.sum())[0]
    return expr


def get_droplet_expr(
    path_rna: Union[AnyStr, pd.Series, Sequence], plot: bool = False, x_arange: Tuple[int, int] = (0, 2000, 10)
) -> Union[pd.DataFrame, Dict[str, np.ndarray]]:
    if not isinstance(path_rna, (str, Path)):
        # TODO: add parallel implementation
        return pd.Series(path_rna).map(get_droplet_expr).apply(pd.Series)
    expr_dict = dict()

    path_raw = Path(path_rna).parent / "raw_feature_bc_matrix.h5"
    try:
        ad_all = read_10x(path_raw)
    except OSError:
        logging.warning(f"\nCouldn't load {path_rna}. Skipping")
        return {"cell": np.nan, "noncell": np.nan, "empty": np.nan}
    df_filt = pd.read_csv(Path(path_rna) / "barcodes.tsv.gz", index_col=0).index
    ad_cel = ad_all[ad_all.obs.index.isin(df_filt)]
    ad_non = ad_all[~ad_all.obs.index.isin(df_filt)]
    ad_emp = ad_non[ad_non.X.sum(axis=1) < ad_cel.X.sum(axis=1).min() / 10]
    expr_dict["cell"] = _pseudo_bulk_vector(ad_cel)
    expr_dict["noncell"] = _pseudo_bulk_vector(ad_non)
    expr_dict["empty"] = _pseudo_bulk_vector(ad_emp)
    if plot:
        for k, _ad in {"cell": ad_cel, "noncell": ad_non, "empty": ad_emp}.items():
            plt.hist(_ad.X.sum(axis=1), bins=np.arange(*x_arange), log=True, alpha=0.5, label=k)
        plt.title(f"droplet hist {str(path_rna).split('/')[-3]}")
        plt.xlabel("UMIs")
        plt.ylabel("droplets")
        plt.xlim(*x_arange[:2])
        plt.legend()
        plt.show()
    return expr_dict


def _subset_expr_by_meta(df_expr, dff):
    return df_expr[df_expr.index.isin(dff["lane"])]


def _group_expr_by_meta(df_expr, df, col):
    return {val: _subset_expr_by_meta(df_expr, dff) for val, dff in df.groupby(col)}


def droplet_expr_corr(
    expr: Union[pd.Series, pd.DataFrame],
    df=None,
    groupby=None,
    cluster=False,
    title="pseudobulk_similarity",
    plot=True,
    n_pcs=1,
):
    """

    Args:
        expr: series of numpy arrays, where each array is the pseudobulk
            expression profile of a sample
        df: metadata dataframe by which to subset the expression series
        groupby:
        cluster:
        title:
        n_pcs:
    Returns:

    """
    if groupby:
        for name, dff_expr in _group_expr_by_meta(expr, df, groupby).items():
            _title = f"{title}-{name}" if title else name
            droplet_expr_corr(dff_expr, df, None, cluster, title=_title)
        return
    if isinstance(expr, pd.DataFrame):
        for col in expr:
            _title = f"{title}-{col}"
            droplet_expr_corr(expr[col], df, groupby, cluster, title=_title)
        return
    elif df is not None:
        expr = _subset_expr_by_meta(expr, df)
    if len(expr) < 2:
        logging.warning(f"{title} ({len(expr)} samples): too few samples for comparison")
        return
    if plot:
        if cluster:
            figsize = (4 + 0.3 * len(expr), 4 + 0.25 * len(expr))
            kwargs = {"figsize": figsize}
        else:
            figsize = (0.25 * len(expr), 0.2 * len(expr))
            _, ax = plt.subplots(figsize=figsize)
            kwargs = {"ax": ax}
            kwargs["ax"].set_title(title)
        corr, plot = pairwise_metric_heat(expr, lambda x, y: pearsonr(x, y)[0], cluster=cluster, **kwargs)
        if cluster:
            plt.title(title)
        if df is not None:
            pass
        output = [corr, plot]
    else:
        output = [
            pairwise_metric(expr, lambda x, y: pearsonr(x, y)[0]),
        ]
