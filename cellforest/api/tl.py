from typing import Optional

from anndata import AnnData
import scanpy as sc


def umap(
    ad: AnnData,
    key_added: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    **kwargs,
):
    sc.tl.umap(ad, neighbors_key=neighbors_key, **kwargs)
    ad.obsm[key_added] = ad.obsm["X_umap"]


def tsne(
    ad: AnnData,
    key_added: Optional[str] = None,
    neighbors_key: Optional[str] = None,
    **kwargs,
):
    sc.tl.tsne(ad, neighbors_key=neighbors_key, **kwargs)
    ad.obsm[key_added] = ad.obsm["X_tsne"]


def reset_colors(ad: AnnData):
    keys = [k for k in ad.uns_keys() if "_colors" in k]
    [ad.uns.pop(k) for k in keys]
