import scanpy as sc
from cellforest.utils.scanpy.cell import reduce


def show_harmony(ad, show_orig=True, color="sample"):
    ad = ad.copy()
    if show_orig:
        sc.pl.umap(ad, color=color)
    sc.external.pp.harmony_integrate(ad, key="sample")
    ad = reduce(ad, neighbors_kwargs={"use_rep": "X_pca_harmony"})
    sc.pl.umap(ad, color=color)
    return ad


def show_bbknn(ad, show_orig=True, color="sample"):
    ad = ad.copy()
    if show_orig:
        sc.pl.umap(ad, color=color)
    ad = sc.external.pp.bbknn(ad, batch_key="sample", copy=True)
    sc.tl.umap(ad)
    sc.pl.umap(ad, color=color)
    return ad


def show_batch_corr(ad, show_orig=True, color="sample"):
    if show_orig:
        sc.pl.umap(ad, color=color)
    print("harmony")
    show_harmony(ad, show_orig=False, color=color)
    print("bbknn")
    show_bbknn(ad, show_orig=False, color=color)
