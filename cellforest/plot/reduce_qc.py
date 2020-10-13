from dataforest.plot import plot_py, plot_r, requires

import matplotlib.pyplot as plt

from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_r_script_logged


@requires("reduce")
@plot_r
def plot_pca_elbow_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_elbow_curv")


@requires("reduce")
@plot_r
def plot_pca_loadings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_loadings_scat")


@requires("reduce")
@plot_r
def plot_pca_embeddings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_embeddings_scat")


@requires("reduce")
@plot_py
def plot_umap_embeddings_scat(branch: "CellBranch", **kwargs):
    ax = kwargs.pop("ax", plt.gca())
    if "stratify" in kwargs:
        for (name, grp) in branch.meta[["sample_id", "UMAP_1", "UMAP_2"]].groupby("sample_id"):
            ax.scatter(grp["UMAP_1"], grp["UMAP_2"], label=name, **kwargs)
    else:
        ax.scatter(branch.meta["UMAP_1"], branch.meta["UMAP_2"], **kwargs)
    ax.legend()


@requires("reduce")
@plot_r
def plot_umap_embeddings_scat_r(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_umap_embeddings_scat")


@requires("reduce")
@plot_r
def plot_cell_cycle_scoring_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_cell_cycle_scoring_scat")
