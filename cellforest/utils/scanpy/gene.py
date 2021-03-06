from typing import List

from anndata import AnnData
import numpy as np
import scanpy as sc

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
