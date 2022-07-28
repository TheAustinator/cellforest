from typing import Optional

import numpy as np
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


def undersample(
    ad: AnnData,
    obs_class: str,
    max_class: Optional[int] = None,
    max_tot: Optional[int] = None,
) -> AnnData:
    """
    Undersample classes in an anndata either by a total max or a per class max.
    """
    if max_tot:
        raise NotImplementedError()
    all_inds = list()
    slicer = np.array(len(ad) * [False])
    max_class = max_class if max_class else ad.obs[obs_class].value_counts().min()
    for class_ in ad.obs[obs_class].unique():
        inds = ad.obs.reset_index()[ad.obs.reset_index()[obs_class] == class_].index
        inds = np.random.choice(inds, min(len(inds), max_class), replace=False).tolist()
        all_inds += inds
    slicer[all_inds] = True
    return ad[slicer]
