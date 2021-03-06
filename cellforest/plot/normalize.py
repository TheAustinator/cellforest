from dataforest.plot import plot_py, plot_r
import matplotlib.pyplot as plt

from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_r_script_logged


@plot_py(requires="normalize")
def plot_perc_mito_per_cell_hist(branch: "CellBranch", **kwargs):
    _hist_col(branch, "percent.mito", **kwargs)


@plot_py(requires="normalize")
def plot_perc_ribo_per_cell_hist(branch: "CellBranch", **kwargs):
    _hist_col(branch, "percent.ribo", **kwargs)


@plot_py(requires="normalize")
def plot_perc_hsp_per_cell_hist(branch: "CellBranch", **kwargs):
    _hist_col(branch, "percent.hsp", **kwargs)


@plot_py(requires="normalize")
def plot_umis_vs_perc_mito_scat(branch: "CellBranch", **kwargs):
    _scatter_umi_vs_col(branch, "percent.mito", **kwargs)


@plot_py(requires="normalize")
def plot_umis_vs_perc_ribo_scat(branch: "CellBranch", **kwargs):
    _scatter_umi_vs_col(branch, "percent.ribo", **kwargs)


@plot_py(requires="normalize")
def plot_umis_vs_perc_hsp_scat(branch: "CellBranch", **kwargs):
    _scatter_umi_vs_col(branch, "percent.hsp", **kwargs)


@plot_r(requires="normalize")
def plot_perc_mito_per_cell_vln_r(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_perc_mito_per_cell_vln_r")


@plot_r(requires="normalize")
def plot_umis_vs_perc_mito_scat_r(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_umis_vs_perc_mito_scat")


@plot_r(requires="normalize")
def plot_perc_ribo_per_cell_vln_r(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_perc_ribo_per_cell_vln_r")


@plot_r(requires="normalize")
def plot_umis_vs_perc_ribo_scat_r(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_umis_vs_perc_ribo_scat")


@plot_r(requires="normalize")
def plot_perc_hsp_per_cell_vln_r(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_perc_hsp_per_cell_vln_r")


def _hist_col(branch: "CellBranch", col, **kwargs):
    ax = kwargs.pop("ax", plt.gca())
    ax.hist(branch.meta[col], **kwargs)
    ax.set_xlabel(col)


def _scatter_umi_vs_col(branch: "CellBranch", col, **kwargs):
    ax = kwargs.pop("ax", plt.gca())
    y_col = "nCount_RNA"
    ax.scatter(branch.meta[col], branch.meta[y_col], **kwargs)
    ax.set_xlabel(col)
    ax.set_ylabel(y_col)
