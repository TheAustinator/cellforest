import math
from typing import Optional, Literal, Union, List

from dataforest.plot import plot_py, plot_r
from dataforest.utils import listify
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
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


@plot_py(forbid=["facet", "stratify"])
def plot_hist_features(
    branch: "CellBranch",
    features: list,
    ncol: int = 3,
    ax_size=5,
    sep: Literal["stratify", "facet"] = "facet",
    **kwargs,
):
    # TODO: add facet_vars to plotting decorator to allow facet by cols rather than rows, and use that for grid
    rna_features = list(set(features).difference(branch.meta.columns))
    if rna_features:
        branch = branch.copy()
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
    fig.set_size_inches(1.1 * ncol * ax_size, nrow * ax_size)  # 1.1 for color bar
    for i, feat in enumerate(features):
        ax = ax_list[i]
        ax.set_title(feat)
        _plot_feature_hist(branch.meta[feat], ax=ax, **kwargs)


@plot_py
def plot_expressing_cells_bar(
    branch: "CellBranch", genes: Union[str, List[str]], x: Optional[str] = None, hue: Optional[str] = None, **kwargs
):
    df = _expressing_cells_df(branch, genes, x, hue)
    return barplot(x="x", y="frac_expr", hue="hue", data=df, **kwargs)


@plot_py
def plot_expressing_cells_bar_facet_gene(
    branch: "CellBranch",
    genes: Union[str, List[str]],
    x: Optional[str] = None,
    hue: Optional[str] = None,
    ncol: int = 4,
    ax_size: int = 5,
    **kwargs,
):

    ncol = min(len(genes), ncol)
    nrow = math.ceil(len(genes) / ncol)
    fig, ax_arr = plt.subplots(nrow, ncol)
    if not isinstance(ax_arr, np.ndarray):
        ax_arr = np.array([[ax_arr]])
    if ax_arr.ndim < 2:
        np.expand_dims(ax_arr, axis=0)
    ax_list = ax_arr.flatten()
    kwargs.pop("ax", None)
    fig.set_size_inches(ncol * ax_size, nrow * ax_size)
    # TODO: will be much faster if
    df = _expressing_cells_df(branch, genes, x, hue)
    for i, (gene, sub_df) in enumerate(df.groupby("gene")):
        ax = ax_list[i]
        ax.set_title(gene)
        return barplot(x="x", y="frac_expr", hue="hue", data=sub_df, **kwargs)


@plot_r
def plot_highest_exprs_dens(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_highest_exprs_dens")


@plot_r
def plot_umis_per_barcode_rank_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, branch.current_process)


def _plot_feature_hist(feature_arr, ax, **kwargs):
    return ax.hist(feature_arr, **kwargs)


def _expressing_cells_df(
    branch: "CellBranch", genes: Union[str, List[str]], x: Optional[str] = None, hue: Optional[str] = None,
):
    def _calc_frac_expressing(branch_, genes_):
        expr_arr = np.array((branch_.rna[:, genes_] > 0).sum(axis=0))[0] / len(branch_.rna)
        return dict(zip(genes_, expr_arr))

    def _get_subset_expr_fracs(name_, br_sub_):
        subset_dicts = list()
        name_ = listify(name_)
        d = dict(zip(grp_cols, name_))
        expr_dict = _calc_frac_expressing(br_sub_, genes)
        for gene, frac_expr in expr_dict.items():
            sub_d = {**d, "gene": gene, "frac_expr": frac_expr}
            subset_dicts.append(sub_d)
        return subset_dicts

    x = listify(x)
    hue = listify(hue)
    genes = listify(genes)
    grp_cols = list(x + hue)
    if "gene" in grp_cols:
        grp_cols.remove("gene")
    # TODO: extremely slow -- parallelize with swifter or pandarallel
    subset_dict_list = [_get_subset_expr_fracs(name, br_sub) for name, br_sub in branch.groupby(by=grp_cols)]
    subset_dicts = [d for l in subset_dict_list for d in l]
    df = pd.DataFrame(subset_dicts)
    df["x"] = df[x].agg(" ".join, axis=1)
    df["hue"] = df[hue].agg(" ".join, axis=1)
    return df
