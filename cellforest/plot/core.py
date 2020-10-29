from typing import Optional

from dataforest.plot import plot_py, plot_r
from seaborn import barplot, violinplot

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


@plot_py
def plot_meta_vln(branch: "CellBranch", **kwargs):
    """
    Args:
        branch:
        **kwargs: for `seaborn.violinplot`
    """
    return violinplot(data=branch.meta, **kwargs)


@plot_py
# TODO: make `plot_py` into class so that you can specify `pass_stratify` to keep stratify
# TODO: if lists for x and hue, join columns as strs to make compatible with seaborn
def plot_frac_cells_recovered_bar(branch: "CellBranch", x: Optional[str] = None, hue: Optional[str] = None, **kwargs):
    """
    Args:
        branch:
        x: see `seaborn.barplot`
        hue: see `seaborn.barplot`
        **kwargs: passed to seaborn
    """
    lane_meta = branch.meta.drop_duplicates()

    def _get_cells_loaded(row):
        qry = " and ".join([f"{k} == '{v}'" for k, v in dict(row).items()])
        meta_sub = lane_meta.query(qry)
        cells_loaded = meta_sub["cells_loaded"].sum()
        return cells_loaded

    grp_bar = [] if not x else [x]
    grp_clr = [] if not hue else [hue]
    grp_key = grp_bar + grp_clr
    col = branch.meta.columns[0]
    df = branch.meta.groupby(grp_key).aggregate(len)
    df = df[col].reset_index().rename(columns={col: "cells_recovered"})
    df["cells_loaded"] = df[grp_key].apply(_get_cells_loaded, axis=1)
    df["frac_recovered"] = df["cells_recovered"] / df["cells_loaded"]
    return barplot(x=x, y="frac_recovered", hue=hue, data=df, **kwargs)


@plot_r
def plot_highest_exprs_dens(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_highest_exprs_dens")


@plot_r
def plot_umis_per_barcode_rank_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, branch.current_process)
