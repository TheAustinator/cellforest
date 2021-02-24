import matplotlib.pyplot as plt
import numpy as np
import scanpy as sc

LINE = {"color": "r", "linewidth": 1}


def volcano(
    df,
    x="logfc",
    y="-logp",
    c="mean_expr_clip",
    group=None,
    pval_thresh=None,
    logfc_thresh=None,
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
        print(name)
        plt.show()


def plot_gene_set(ad, df, gene_set, alpha=1, s=10, groupby="cell_type", volcano_groups=None, **kwargs):
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
    volcano_groups = volcano_groups if volcano_groups is not None else ad.obs[groupby].unique()
    volcano_groups = volcano_groups if not isinstance(volcano_groups, (str, int)) else (volcano_groups,)
    for group in volcano_groups:
        volcano(
            df[df["gene"].isin(gene_set)], group=group, alpha=alpha, s=s, **kwargs, logfc_thresh=0.25, pval_thresh=0.01
        )
        df_gs = df[(df["gene"].isin(gene_set)) & (df["group"] == group)]
        df_gs = df_gs.sort_values(["mean_expr"])
        sc.pl.heatmap(ad, var_names=df_gs["gene"], groupby=groupby, standard_scale="var")
