from dataforest.plot import plot_py, plot_r

from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_r_script_logged


@plot_py
def plot_genes_per_cell_hist(branch: "CellBranch", **kwargs):
    branch.rna.hist("nonzero", axis=0, **kwargs)


@plot_py
def plot_umis_per_cell_hist(branch: "CellBranch", **kwargs):
    branch.rna.hist("sum", axis=0, **kwargs)


@plot_py
def plot_umis_vs_genes_scat(branch: "CellBranch", **kwargs):
    branch.rna.scatter(agg_x="nonzero", agg_y="sum", axis=0, **kwargs)


@plot_r
def plot_highest_exprs_dens(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_highest_exprs_dens")


@plot_r
def plot_umis_per_barcode_rank_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, branch.current_process)
