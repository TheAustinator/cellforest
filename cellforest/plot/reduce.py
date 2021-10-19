import math
from collections import Iterable, Sized
from typing import List

from dataforest.plot import plot_py, plot_r
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_r_script_logged


@plot_py(requires="reduce")
def plot_umap_embeddings_scat(branch: "CellBranch", **kwargs):
    _plot_umap_embeddings_scat(branch, **kwargs)


@plot_py(requires="reduce", forbid=["facet", "stratify"])
def plot_umap_features(
    branch: "CellBranch", features: List, ncol: int = 3, ax_size=5, **kwargs
):
    # TODO: add facet_vars to plotting decorator to allow facet by cols rather than rows, and use that for grid
    rna_features = list(set(features).difference(branch.meta.columns))
    if rna_features:
        branch = branch.copy()
        missing = set(rna_features).difference(branch.rna.genes.values)
        if missing:
            raise ValueError(f"Genes not found: {missing}")
        branch.meta[rna_features] = branch.rna[:, rna_features].todense()
    ncol = min(len(features), ncol)
    nrow = math.ceil(len(features) / ncol)
    fig, ax_arr = plt.subplots(nrow, ncol)
    if not isinstance(ax_arr, np.ndarray):
        ax_arr = np.array([[ax_arr]])
    if ax_arr.ndim < 2:
        np.expand_dims(ax_arr, axis=0)
    ax_list = ax_arr.flatten()
    kwargs.pop("ax", None)
    fig.set_size_inches(1.1 * ncol * ax_size, nrow * ax_size)  # 1.1 for key bar
    for i, feat in enumerate(features):
        ax = ax_list[i]
        ax.set_title(feat)
        scat = _plot_umap_embeddings_scat(branch, c=branch.meta[feat], ax=ax, **kwargs)
        fig.colorbar(scat, ax=ax, alpha=1, pad=0)


@plot_r(requires="reduce")
def plot_pca_elbow_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_elbow_curv")


@plot_r(requires="reduce")
def plot_pca_loadings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_loadings_scat")


@plot_r(requires="reduce")
def plot_pca_embeddings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_embeddings_scat")


@plot_r(requires="reduce")
def plot_umap_embeddings_scat_r(
    branch: "CellBranch", r_script: str, args: list, **kwargs
):
    run_r_script_logged(branch, r_script, args, "plot_umap_embeddings_scat")


@plot_r(requires="reduce")
def plot_cell_cycle_scoring_scat(
    branch: "CellBranch", r_script: str, args: list, **kwargs
):
    run_r_script_logged(branch, r_script, args, "plot_cell_cycle_scoring_scat")


def _plot_umap_embeddings_scat(branch: "CellBranch", **kwargs):
    ax = kwargs.pop("ax", plt.gca())
    if "stratify" in kwargs:
        for (name, grp) in branch.meta[["sample_id", "UMAP_1", "UMAP_2"]].groupby(
            "sample_id"
        ):
            scat = ax.scatter(grp["UMAP_1"], grp["UMAP_2"], label=name, **kwargs)
    else:
        scat = ax.scatter(branch.meta["UMAP_1"], branch.meta["UMAP_2"], **kwargs)
    ax.legend()
    return scat
