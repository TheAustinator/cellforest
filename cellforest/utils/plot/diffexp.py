import matplotlib.pyplot as plt
import numpy as np

LINE = {"color": "r", "linewidth": 1}


def volcano(
    df,
    x="logfc",
    y="-logp",
    c="mean_expr_clip",
    group=None,
    pval_thresh=0.5,
    logfc_thresh=0.25,
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
        logp_thresh = -np.log10(pval_thresh)
        plt.axvline(logfc_thresh, **LINE)
        plt.axvline(-logfc_thresh, **LINE)
        plt.axhline(logp_thresh, **LINE)
        if xlim:
            plt.xlim(*xlim)
        if ylim:
            plt.ylim(*ylim)
        plt.title(name)
        plt.ylabel("-logp")
        plt.xlabel("logfc")
        print(name)
        plt.show()
