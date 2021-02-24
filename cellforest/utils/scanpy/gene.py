from typing import List

from anndata import AnnData
import numpy as np
import scanpy as sc
from sklearn.preprocessing import normalize, scale

from cellforest.utils.scanpy.generic import _generic_preprocess


def drop_gene_prefix(ad: AnnData, prefix_list: List[str]):
    var_names = ad.var_names[~ad.var_names.str.startswith(prefix_list)]
    ad = ad[:, var_names]
    return ad


def agg_gene_prefix(ad: AnnData, prefix_list: List[str], new_obs_colname: str, drop: bool = True):
    var_names = ad.var_names[ad.var_names.str.startswith(prefix_list)]
    ad.obs[new_obs_colname] = ad[:, var_names].X.sum(axis=1)
    if drop:
        ad = drop_gene_prefix(ad, prefix_list)
    return ad


def add_prefix_fracs(ad, prefix_list):
    for prefix in prefix_list:
        ad.obs[prefix] = ad[:, ad.var.index.str.startswith(prefix)].X.mean(axis=1) / ad.X.mean(axis=1)
    return ad


def process_transpose(
    ad: AnnData, min_cells: int = 300, min_genes: int = 200, max_genes: int = 2500, max_pct_mito: int = 30
):
    ad = ad.copy()
    ad = _generic_preprocess(ad, min_cells, min_genes, max_genes, max_pct_mito)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, batch_key="sample")
    ad = ad.transpose()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.pca(ad, n_comps=50)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    return ad


def get_genes_clustered(adt, markers):
    sc.tl.umap(adt, min_dist=0.00001, spread=1)
    sc.tl.leiden(adt, resolution=5)
    sc.pl.umap(adt, color="leiden")
    marker_clusters = adt.obs[adt.obs_names.isin(markers)]["leiden"]
    print(marker_clusters)
    marker_umap_coords = np.array(adt[adt.obs_names.isin(markers)].obsm["X_umap"])
    print(marker_umap_coords)
    clusters = marker_clusters.unique()
    genes = list(adt.obs_names[adt.obs["leiden"].isin(clusters)])
    return genes


def filter_markers(df, logfc=0, pval_adj=1, mean_expr=0, frac_expr=0, filter_prefixes=(), group=None):
    df = df[
        (df["logfc"].abs() > logfc)
        & (df["pval_adj"] < pval_adj)
        & (df["mean_expr"] > mean_expr)
        & (df["frac_expr"] > frac_expr)
    ]
    if group:
        df = df[df["group"] == group]
    df = df[~df["gene"].str.startswith(filter_prefixes)]
    return df


def rank_markers(df):
    return df["logfc"] * -np.log10(df["pval_adj"])


def get_feature(ad, f, std_scale=False, norm=False):
    if not isinstance(f, str):
        return np.vstack(list(map(lambda _f: get_feature(ad, _f), f))).T
    if f in ad.obs.columns:
        arr = ad.obs[f]
    else:
        arr = ad[:, f].X.toarray().T[0]
    if std_scale:
        arr = scale(arr.T).T
    if norm:
        arr = normalize(arr.T).T
    return arr
