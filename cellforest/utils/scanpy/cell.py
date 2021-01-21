from functools import wraps

from anndata import AnnData
from dataforest.utils.decorators import sub_kwargs
import numpy as np
import pandas as pd
import scanpy as sc

from cellforest.utils.scanpy.generic import _generic_preprocess


# maybe increase `min_cells
def preprocess(ad: AnnData, min_cells: int = 10, min_genes: int = 200, max_genes: int = 2500, max_pct_mito: int = 30):
    ad = _generic_preprocess(ad, min_cells, min_genes, max_genes, max_pct_mito)
    ad.raw = ad
    sc.pp.log1p(ad)
    return ad


@sub_kwargs("hvg_kwargs", "pca_kwargs", "neighbors_kwargs", "umap_kwargs")
def reduce(ad: AnnData, n_comps: int = 15, hvg_kwargs=None, pca_kwargs=None, neighbors_kwargs=None, umap_kwargs=None):
    sc.pp.highly_variable_genes(ad, batch_key="sample", **hvg_kwargs)
    sc.pp.pca(ad, n_comps=n_comps, **pca_kwargs)
    sc.pp.neighbors(ad, **neighbors_kwargs)
    sc.tl.umap(ad, **umap_kwargs)
    return ad


def cluster(ad: AnnData, resolution: float = 0.1, n_genes: int = 20, fontsize: int = 8, ax=None):
    sc.tl.leiden(ad, resolution=resolution)
    markers(ad, "leiden", n_genes, fontsize, ax)
    return ad


def markers(ad: AnnData, key="leiden", n_genes: int = 20, fontsize: int = 8, ax=None):
    sc.pl.umap(ad, color=key)
    sc.tl.rank_genes_groups(ad, groupby=key)
    sc.pl.rank_genes_groups(ad, n_genes=n_genes, fontsize=fontsize, ax=ax)
    return ad


def process(ad: AnnData):
    ad = preprocess(ad)
    ad = reduce(ad)
    return ad


def get_markers_df(ad: AnnData, mean_expr_clip: float = 0.1, group=None):
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
        try:
            mean_expr = ad[ad.obs[group_var] == _group, _dff["gene"]].X.mean(axis=0)
        except Exception as e:
            print(ad[ad.obs[group_var] == _group, _dff["gene"]].X.shape)
            raise e
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
    return df
