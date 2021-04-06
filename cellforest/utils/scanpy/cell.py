from functools import wraps

from anndata import AnnData
from dataforest.utils.decorators import sub_kwargs
import numpy as np
import pandas as pd
import scanpy as sc

from cellforest.utils.scanpy.generic import _generic_preprocess


# maybe increase `min_cells
def preprocess(
    ad: AnnData, min_cells: int = None, min_genes: int = None, max_genes: int = None, max_pct_mito: int = None
):
    # if "raw" not in ad.layers:
    ad.layers["raw"] = ad.X.copy()
    ad.raw = ad.copy()
    ad = _generic_preprocess(ad, min_cells, min_genes, max_genes, max_pct_mito)
    sc.pp.log1p(ad)
    return ad


@sub_kwargs("hvg_kwargs", "pca_kwargs", "neighbors_kwargs", "umap_kwargs", "tsne_kwargs")
def reduce(
    ad: AnnData,
    n_comps: int = 15,
    hvg_kwargs=None,
    pca_kwargs=None,
    neighbors_kwargs=None,
    umap_kwargs=None,
    tsne_kwargs=None,
    umap=True,
    tsne=False,
):
    sc.pp.highly_variable_genes(ad, batch_key="sample", **hvg_kwargs)
    sc.pp.pca(ad, n_comps=n_comps, **pca_kwargs)
    sc.pp.neighbors(ad, **neighbors_kwargs)
    if umap:
        sc.tl.umap(ad, **umap_kwargs)
    if tsne:
        sc.tl.tsne(ad, **tsne_kwargs)
    return ad


def cluster(ad: AnnData, resolution: float = 0.1, **markers_kwargs):
    sc.tl.leiden(ad, resolution=resolution)
    markers(ad, "leiden", **markers_kwargs)
    return ad


def markers(
    ad: AnnData, key="leiden", n_genes: int = 20, filter_var_prefix=(), obs_facet=None, fontsize: int = 8, ax=None,
):
    if filter_var_prefix:
        ad = ad[:, ~ad.var_names.str.startswith(filter_var_prefix)]
        del ad.raw
    if obs_facet:
        for val, _obs in ad.obs.groupby(obs_facet):
            print(val)
            _ad = ad[ad.obs_names.isin(_obs.index)]
            sc.pl.umap(_ad, color=key)
    else:
        sc.pl.umap(ad, color=key)
    sc.tl.rank_genes_groups(ad, groupby=key)
    sc.pl.rank_genes_groups(ad, n_genes=n_genes, fontsize=fontsize, ax=ax)
    get_markers_df(ad)
    return ad


def process(ad: AnnData):
    ad = preprocess(ad)
    ad = reduce(ad)
    return ad


def get_markers_df(ad: AnnData, mean_expr_clip: float = 0.1, group=None, uns_key="markers"):
    """

    Args:
        ad:
        mean_expr_clip:
        group:

    Returns:

    """

    def get_col(col):
        return np.array(ad.uns["rank_genes_groups"][col].tolist())

    cols = ["names", "pvals_adj", "logfoldchanges"]
    col_dict = {col: get_col(col).flatten() for col in cols}
    labels = ad.uns["rank_genes_groups"]["names"].dtype.names
    size = len(get_col("names"))
    groups = np.tile(labels, (size, 1)).flatten()
    col_dict["group"] = groups
    df = pd.DataFrame(col_dict)
    df.rename(columns={"names": "gene", "pvals_adj": "pval_adj", "logfoldchanges": "logfc"}, inplace=True)
    group_var = ad.uns["rank_genes_groups"]["params"]["groupby"]

    def _add_group(_group, _dff):
        selector = ad.obs[group_var] == _group
        if selector.sum() == 0:
            raise ValueError(
                f"specified `group`: {group}   not in column :{group_var}  . Choose from {ad.obs[group_var].unique()}"
            )
        mean_expr = ad[ad.obs[group_var] == _group, _dff["gene"]].X.mean(axis=0)
        frac_expr = (ad[ad.obs[group_var] == _group, _dff["gene"]].X > 0).mean(axis=0)
        df.loc[_dff.index, "mean_expr"] = mean_expr.T
        df.loc[_dff.index, "frac_expr"] = frac_expr.T

    if not group:
        for group, dff in df.groupby("group"):
            _add_group(group, dff)
    else:
        _add_group(group, df)
    df["mean_expr_clip"] = df["mean_expr"].apply(lambda x: np.clip(x, 0, mean_expr_clip))
    df["-logp"] = -df["pval_adj"].apply(np.log10)
    floor = df["-logp"][df["-logp"].abs() != np.inf].min()
    df["-logp"] = df["-logp"].replace({-np.inf: floor})
    if uns_key:
        ad.uns[uns_key] = df
    return df
