import matplotlib
import matplotlib.pyplot as plt

from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_process_r_script
from cellforest.plot.qc_plot_wrappers import qc_plot_py, qc_plot_r

matplotlib.use("Agg")  # plotting in the backend


@qc_plot_py
def plot_genes_per_cell_hist(branch: "CellBranch", **kwargs):
    branch.rna.hist("nonzero", axis=0, **kwargs)


@qc_plot_py
def plot_umis_per_cell_hist(branch: "CellBranch", **kwargs):
    branch.rna.hist("sum", axis=0, **kwargs)


@qc_plot_py
def plot_umis_vs_genes_hist(branch: "CellBranch", **kwargs):
    branch.rna.scatter(agg_x="nonzero", agg_y="sum", axis=0, **kwargs)


@qc_plot_r
def plot_perc_mito_per_cell_vln(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_process_r_script(branch, r_script, args, branch.current_process)


@qc_plot_r
def plot_highest_expression_genes_dens(branch: "CellBranch", r_script: str, args: list, **kwargs):
    pass
    # R script

