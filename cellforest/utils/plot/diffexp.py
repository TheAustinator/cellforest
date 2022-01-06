import logging
from itertools import product

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from cellforest.utils.scanpy.cell import markers
from cellforest.utils.scanpy.group import groupby_dict

LINE = {"key": "r", "linewidth": 1}


def volcano(
    df,
    x="logfc",
    y="-logp",
    c="mean_expr_clip",
    group=None,
    pval_thresh=None,
    logfc_thresh=None,
    volcano_thresh=None,
    mean_expr_thresh=None,
    frac_expr_thresh=None,
    xlim=None,
    ylim=None,
    s=2,
    alpha=0.2,
    **kwargs,
):
    if group:
        df = df[df["group"] == group]
    for name, dff in df.groupby("group"):
        if isinstance(c, (list, set, tuple)):
            color = dff["gene"].isin(c)
        else:
            color = dff[c]
        plt.scatter(dff[x], dff[y], c=color, s=s, alpha=alpha, **kwargs)
        if pval_thresh is not None:
            logp_thresh = -np.log10(pval_thresh)
            plt.axhline(logp_thresh, **LINE)
        if logfc_thresh is not None:
            plt.axvline(logfc_thresh, **LINE)
            plt.axvline(-logfc_thresh, **LINE)
        if xlim is not None:
            plt.xlim(*xlim)
        if ylim is not None:
            plt.ylim(*ylim)
        plt.title(name)
        plt.ylabel("-logp")
        plt.xlabel("logfc")
        plt.colorbar()
        if volcano_thresh:
            ax = plt.gca()
            xlim = ax.get_xlim()
            ylim = ax.get_ylim()
            range_ = (0.001, xlim[1]) if volcano_thresh > 0 else (xlim[0], -0.001)
            x = np.linspace(*range_, 1000)
            y = volcano_thresh / x
            ax.plot(x, y, c="y")
            ax.set_xlim(xlim)
            ax.set_ylim(ylim)
            return ax
        print(name)
        plt.show()


def plot_gene_set(
    ad, df, gene_set, alpha=1, s=10, groupby="cell_type", volcano_groups=None, **kwargs
):
    """

    Args:
        ad:
        df:
        gene_set:
        alpha:
        s:
        groupby:
        **kwargs:

    Returns:

    """
    volcano_groups = (
        volcano_groups if volcano_groups is not None else ad.obs[groupby].unique()
    )
    volcano_groups = (
        volcano_groups
        if not isinstance(volcano_groups, (str, int))
        else (volcano_groups,)
    )
    for group in volcano_groups:
        volcano(
            df[df["gene"].isin(gene_set)],
            group=group,
            alpha=alpha,
            s=s,
            **kwargs,
            logfc_thresh=0.25,
            pval_thresh=0.01,
        )
        df_gs = df[(df["gene"].isin(gene_set)) & (df["group"] == group)]
        df_gs = df_gs.sort_values(["mean_expr"])
        sc.pl.heatmap(
            ad, var_names=df_gs["gene"], groupby=groupby, standard_scale="var"
        )


def pairwise_de_corr(ad, grp_obs, markers_df, plot=True):
    de_corr = pd.DataFrame()
    de_size = pd.DataFrame()
    ad_grp = groupby_dict(ad, grp_obs)
    combos = list(product(ad_grp.items(), ad_grp.items()))
    markers_grps = markers_df["group"].unique()
    if len(markers_grps) > 1:
        grp_keep = markers_grps[0]
        logging.warning(f"Multiple DE groups: {markers_grps}, keeping {grp_keep}")
        markers_df = markers_df[markers_df["group"] == grp_keep]
    volc = markers_df.set_index("gene")["volcano"]
    for (name_1, ad_1), (name_2, ad_2) in combos:
        _ad = ad_1.concatenate(ad_2, join="outer")
        if len(_ad.obs[grp_obs].unique()) < 2:
            de_corr.loc[name_1, name_2] = 1
            de_size.loc[name_1, name_2] = 1
            continue
        markers(_ad, grp_obs, plot=plot)
        mark = _ad.uns["markers"].set_index("gene")
        mark = mark[mark["group"] == name_1]

        volc_df = pd.DataFrame({"volc_main": volc, "volc_subset": mark["volcano"]})
        corr = volc_df.corr().to_dict()["volc_main"]["volc_subset"]
        size = volc_df["volc_subset"].abs().sum() / volc_df["volc_subset"].abs().sum()
        de_corr.loc[name_1, name_2] = corr
        de_size.loc[name_1, name_2] = size
    sns.heatmap(de_corr)
    plt.title("de_corr")
    plt.show()
    sns.heatmap(de_size)
    plt.title("de_size")
    plt.show()
    return de_corr, de_size
