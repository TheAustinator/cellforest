from copy import deepcopy
from typing import Optional, Union, Literal, List, Collection

from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
from sklearn.utils import sparsefuncs
from scipy.sparse import issparse, hstack

from cellforest.api.tl import tsne, umap
from cellforest.api.pl import embedding

ASSAYS = ("rna", "fb", "crispr", "custom")
METHODS = ("harmony", "bbknn", "scanorama", "mnn")
ASSAY_KEYS = {
    "rna": "Gene Expression",
    "fb": "Antibody Capture",
    "crispr": "CRISPR Guide Capture",
    "custom": "Custom",
}
ASSAY_PER_FEATURE = {"rna": 0.2732, "fb": 50, "crispr": 1, "custom": 1}
METHOD_TYPES = {
    "harmony": "embedding",
    "scanorama": "embedding",
    "bbknn": "neighbors",
    "combat": "matrix",
}
REDUCERS = {
    "tsne": tsne,
    "umap": umap,
}


def norm_total_per_assay(
    ad: AnnData,
    target_sum: Optional[Union[int, dict]] = None,
    key_added: Optional[str] = "norm",
    assay_keys: Optional[dict] = None,
):
    ad.var_names_make_unique()
    if "raw" in ad.layers:
        ad.X = ad.layers["raw"]
    ad.layers["raw"] = ad.X.copy()
    if len(ad.var["feature_types"].unique()) < 2:
        print("single assay -- using scanpy norm")
        sc.pp.normalize_total(ad, target_sum=target_sum, key_added=key_added)
        ad.layers["norm"] = ad.X.copy()
        return
    print("multiple assays -- using custom norm")
    assay_keys = assay_keys if assay_keys else ASSAY_KEYS
    if not target_sum:
        feature_counts = ad.var["feature_types"].value_counts().to_dict()
        target_sum = {
            k: feature_counts[v] * ASSAY_PER_FEATURE[k]
            for k, v in ASSAY_KEYS.items()
            if v in feature_counts
        }
    elif isinstance(target_sum, (int, float)):
        target_sum = {k: target_sum for k in assay_keys}

    else:
        row_tot = ad.X[0, :].sum()
        if row_tot != round(row_tot):
            if ad.raw:
                ad.X = ad.raw
            else:
                raise ValueError(f".X counts have decimals and are not raw")
    feature_types = ad.var["feature_types"].unique()
    X_assays = {
        k: ad[:, ad.var["feature_types"] == v].X.copy()
        for k, v in assay_keys.items()
        if v in feature_types
    }
    var_assays = {
        k: ad.var[ad.var["feature_types"] == v]
        for k, v in assay_keys.items()
        if v in feature_types
    }
    mat_list = list()
    for assay, X in X_assays.items():
        print(f"norming: {assay}")
        counts = np.asarray(X.sum(axis=1))
        ad.obs[f"{assay}_counts"] = counts
        counts += counts == 0
        counts = counts / target_sum[assay]
        if issparse(X):
            sparsefuncs.inplace_row_scale(X, 1 / counts)
        else:
            np.divide(X, counts[:, None], out=X)
        mat_list.append(X)

    X = hstack(mat_list).tocsr()
    var_names = pd.concat(var_assays.values()).index.tolist()
    ad = ad[:, var_names]
    if key_added:
        ad.layers[key_added] = X.copy()
    ad.X = X


def pca(
    ad: AnnData,
    key_added: Optional[str] = None,
    assays: Union[
        Literal["rna", "fb", "crispr" "custom"],
        List[Literal["rna", "fb", "crispr" "custom"]],
    ] = "rna",
    var_names: Optional[List[str]] = None,
    vars_excl: Optional[List[str]] = None,
    layer: Optional[str] = None,
    **kwargs,
):
    assays = (
        assays if isinstance(assays, (list, set, tuple)) else [assays,]
    )
    assays = sorted(assays)
    key_added = key_added if key_added else f"X_pca_{'_'.join(assays)}"
    feature_types = [ASSAY_KEYS[a] for a in assays]
    ad_assay = ad[:, ad.var["feature_types"].isin(feature_types)]
    if var_names:
        ad_assay = ad_assay[:, var_names]
    if vars_excl:
        ad_assay = ad_assay[: ~ad_assay.var_names.isin(vars_excl)]
    sc.pp.scale(ad_assay, layer=layer)

    sc.pp.pca(ad_assay, **kwargs)
    ad.obsm[key_added] = ad_assay.obsm["X_pca"]


def preprocess_default(ad: AnnData):
    ad.var_names_make_unique()
    if not ad.raw:
        ad.raw = ad
    sc.pp.calculate_qc_metrics(ad, inplace=True)
    norm_total_per_assay(ad)
    sc.pp.log1p(ad)
    ad.layers["log"] = ad.X


def recipe_batch_correct(
    ad: AnnData,
    keys: Union[str, List[str]],
    *,
    assays: Union[
        Literal["rna", "fb", "crispr" "custom"],
        List[Literal["rna", "fb", "crispr" "custom"]],
    ] = "rna",
    methods: Union[
        Literal["harmony", "bbknn", "scanorama", "mnn"],
        List[Literal["harmony", "bbknn", "scanorama", "mnn"]],
    ] = ("harmony", "bbknn", "scanorama"),
    layer: Optional[str] = "log",
    reductions: Union[Literal["tsne", "umap"], List[Literal["tsne", "umap"]]] = "umap",
    force: bool = False,
    plot: bool = True,
    plot_mode: Literal["cf", "sc"] = "sc",
    **plot_kwargs,
):
    if not isinstance(keys, str):
        return [
            recipe_batch_correct(
                ad,
                k,
                assays=assays,
                methods=methods,
                layer=layer,
                reductions=reductions,
                force=force,
                plot=plot,
                plot_mode=plot_mode,
                **plot_kwargs,
            )
            for k in keys
        ]
    assays = (
        assays if isinstance(assays, (list, set, tuple)) else [assays,]
    )
    reductions = (
        reductions if isinstance(reductions, (list, set, tuple)) else [reductions,]
    )
    assay_key = "_".join(sorted(assays))
    assay_batch_key = f"{assay_key}_{keys}"
    pca_key = f"X_pca_{assay_key}"
    harm_key = f"X_harmony_{assay_batch_key}"
    scan_key = f"X_scanorama_{assay_batch_key}"
    bbknn_key = f"bbknn_{assay_batch_key}"
    sc.pp.highly_variable_genes(ad, batch_key=keys, layer=layer)
    pca(ad, key_added=pca_key, assays=assays, layer=layer)
    sc.pp.neighbors(ad, use_rep=pca_key, key_added="orig")
    umap(ad, key_added="X_umap_orig", neighbors_key="orig")
    if "combat" in methods:
        raise NotImplementedError("combat")
    if "mnn" in methods:
        raise NotImplementedError("mnn")
    if "harmony" in methods:
        if harm_key in ad.obsm and not force:
            print(f"{harm_key} already exists. Skipping")
        else:
            sce.pp.harmony_integrate(
                ad, key=keys, basis=pca_key, adjusted_basis=harm_key
            )
    if "scanorama" in methods:
        if scan_key in ad.obsm and not force:
            print(f"{scan_key} already exists. Skipping")
        else:
            _ad = ad[:, ad.var["highly_variable"]].copy()
            sce.pp.scanorama_integrate(
                _ad, key=keys, basis=pca_key, adjusted_basis=scan_key,
            )
            ad.obsm[scan_key] = _ad.obsm[scan_key]
    if "bbknn" in methods:
        if f"connectivities_{bbknn_key}" in ad.obsp and not force:
            print(f"connectivities_{bbknn_key} already exists. Skipping")
        else:
            sc.external.pp.bbknn(ad, batch_key=keys, use_rep=pca_key)
            ad.obsp[f"distances_{bbknn_key}"] = ad.obsp[f"distances"]
            ad.obsp[f"connectivities_{bbknn_key}"] = ad.obsp[f"connectivities"]
            ad.uns[f"neighbors_{bbknn_key}"] = deepcopy(ad.uns["neighbors"])
    for method in methods:
        method_key = f"{method}_{assay_batch_key}"
        neighbors_key = f"neighbors_{method_key}"
        if METHOD_TYPES[method] != "neighbors":
            sc.pp.neighbors(ad, use_rep=f"X_{method_key}", key_added=neighbors_key)
        for redux in reductions:
            reducer = REDUCERS[redux]
            redux_key = f"X_{redux}_{method_key}"
            reducer(ad, key_added=redux_key, neighbors_key=neighbors_key)
            if plot:
                plotter = embedding if plot_mode == "cf" else sc.pl.embedding
                plotter(ad, redux_key, color=keys, **plot_kwargs)
