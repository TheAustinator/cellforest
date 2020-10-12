import matplotlib.pyplot as plt

from cellforest import CellBranch
from cellforest.plot.qc_plot_wrappers import qc_plot_py, requires


@requires("normalize")
@qc_plot_py
def plot_perc_mito_per_cell_hist(branch: "CellBranch", **kwargs):
    _hist_col(branch, "percent.mito", **kwargs)


@requires("normalize")
@qc_plot_py
def plot_perc_ribo_per_cell_hist(branch: "CellBranch", **kwargs):
    _hist_col(branch, "percent.ribo", **kwargs)


@requires("normalize")
@qc_plot_py
def plot_perc_hsp_per_cell_hist(branch: "CellBranch", **kwargs):
    _hist_col(branch, "percent.hsp", **kwargs)


@requires("normalize")
@qc_plot_py
def plot_umis_vs_perc_mito_scat(branch: "CellBranch", **kwargs):
    _scatter_umi_vs_col(branch, "percent.mito", **kwargs)


@requires("normalize")
@qc_plot_py
def plot_umis_vs_perc_ribo_scat(branch: "CellBranch", **kwargs):
    _scatter_umi_vs_col(branch, "percent.ribo", **kwargs)


@requires("normalize")
@qc_plot_py
def plot_umis_vs_perc_mito_scat(branch: "CellBranch", **kwargs):
    _scatter_umi_vs_col(branch, "percent.hsp", **kwargs)


def _hist_col(branch: "CellBranch", col, **kwargs):
    ax = kwargs.pop("ax", plt.gca())
    ax.hist(branch.meta[col], **kwargs)
    ax.set_xlabel(col)


def _scatter_umi_vs_col(branch: "CellBranch", col, **kwargs):
    ax = kwargs.pop("ax", plt.gca())
    ax.scatter(branch.meta[col], branch.meta[(y_col := "nCount_RNA")], **kwargs)
    ax.set_xlabel(col)
    ax.set_ylabel(y_col)
