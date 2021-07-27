from typing import List

from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import vstack, hstack, coo_matrix

from cellforest.utils.scanpy.generic import _generic_preprocess


def drop_gene_prefix(ad: AnnData, prefix_list: List[str]):
    prefixes = "|".join(map(lambda s: f"^{s}", prefix_list))
    var_names = ad.var_names[~ad.var_names.str.contains(prefixes)]
    ad = ad[:, var_names]
    return ad


def agg_gene_prefix(ad: AnnData, prefix_list: List[str], new_obs_colname: str, drop: bool = False):
    prefixes = "|".join(map(lambda s: f"^{s}", prefix_list))
    var_names = ad.var_names[ad.var_names.str.contains(prefixes)]
    ad.obs[new_obs_colname] = ad[:, var_names].X.sum(axis=1)
    if drop:
        ad = drop_gene_prefix(ad, prefix_list)
    return ad


def add_prefix_fracs(ad, prefix_list, raw=True):
    prefix_list = [prefix_list,] if isinstance(prefix_list, (str, int)) else prefix_list
    for prefix in prefix_list:
        ad_tot = ad.raw if raw and ad.raw else ad
        ad_sub = ad[:, ad.var.index.str.startswith(prefix)]
        X_tot = ad_tot.X
        X_sub = ad_sub.X
        ad.obs[prefix] = X_sub.sum(axis=1) / X_tot.sum(axis=1)
    return ad


def process_transpose(
    ad: AnnData, min_cells: int = 10, min_genes: int = 200, max_genes: int = 2500, max_pct_mito: int = 30
):
    ad = ad.copy()
    ad.X = ad.raw.X
    ad = _generic_preprocess(ad, min_cells, min_genes, max_genes, max_pct_mito)
    sc.pp.log1p(ad)
    sc.pp.highly_variable_genes(ad, batch_key="sample")
    ad = ad.transpose()
    sc.pp.pca(ad, n_comps=50)
    sc.pp.neighbors(ad)
    sc.tl.umap(ad)
    return ad


def cluster(adt, min_dist=0.00001, spread=1, resolution=5, plot=True):
    # TODO: change to take in original anndata, transpose, then add these annotations to var
    sc.tl.umap(adt, min_dist=min_dist, spread=spread)
    sc.tl.leiden(adt, resolution=resolution)
    if plot:
        sc.pl.umap(adt, color="leiden")


def query_module(adt, markers):
    markers = [markers,] if isinstance(markers, str) else markers
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


def get_features(ad, f, f_obsm=None, n_obsm=None, ignore_missing: bool = False, std_scale=False, norm=False):
    if f is None:
        return pd.DataFrame()
    ad = ad.copy()
    f = [f] if isinstance(f, str) else f
    if f_obsm:
        ad = add_obsm(ad, f_obsm, n_obsm, copy=True)
        f += ad.obs.columns[ad.obs.columns.str.startswith("|".join(f_obsm))].tolist()
    f_obs = [x for x in f if x in ad.obs.columns]
    f_var = [x for x in f if x in ad.var_names]
    f_mis = set(f).difference(set(f_obs).union(f_var))
    if f_mis and not ignore_missing:
        print(ignore_missing)
        raise ValueError(f"Missing features: {f_mis}. Set `ignore_missing` to ignore.")
    obs = ad.obs[f_obs]
    var = ad[:, f_var].X.toarray() if f_var else []
    var = pd.DataFrame(var, columns=f_var, index=ad.obs.index)
    df = pd.concat([obs, var], axis=1)
    if ignore_missing:
        f = [x for x in f if x in df]
    df = df[f]
    return df


def add_obsm(ad, key_obsm, n_obsm=None, copy=True):
    if isinstance(key_obsm, str):
        key_obsm = [key_obsm]
    for key in key_obsm:
        ad = ad.copy() if copy else ad
        obsm = ad.obsm[key]
        if n_obsm:
            obsm = obsm[:, :n_obsm]
        cols = [key + f"_{i}" for i in range(obsm.shape[1])]
        ad.obs[cols] = obsm
    if copy:
        return ad
