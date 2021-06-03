import scanpy as sc
from cellforest.utils.scanpy.cell import reduce


def harmony(ad, show_orig=True, key="sample", reduce_kwargs=None, **kwargs):
    reduce_kwargs = reduce_kwargs if reduce_kwargs else dict()
    ad = ad.copy()
    if show_orig:
        sc.pl.umap(ad, color=key)
    sc.external.pp.harmony_integrate(ad, key=key, **kwargs)
    ad = reduce(ad, neighbors_kwargs={"use_rep": "X_pca_harmony"}, **reduce_kwargs)
    sc.pl.umap(ad, color=key)
    return ad


def bbknn(ad, show_orig=True, key="sample"):
    ad = ad.copy()
    if show_orig:
        sc.pl.umap(ad, color=key)
    ad = sc.external.pp.bbknn(ad, batch_key=key, copy=True)
    sc.tl.umap(ad)
    sc.pl.umap(ad, color=key)
    return ad


def show_batch_corr(ad, show_orig=True, key="sample"):
    if show_orig:
        sc.pl.umap(ad, color=key)
    print("harmony")
    harmony(ad, show_orig=False, key=key)
    print("bbknn")
    bbknn(ad, show_orig=False, key=key)
