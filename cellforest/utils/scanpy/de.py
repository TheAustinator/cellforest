from frozendict import frozendict
import numpy as np
import pandas as pd
import seaborn as sns


def col_corr(df, col, pivot="donor"):
    df = df[[col, pivot]].pivot(columns=pivot)
    df.columns = df.columns.get_level_values(1)
    return df.corr()


def clustermap_pivot(
    df,
    cols=("logfc", "pval_adj", "-logp", "mean_expr", "mean_expr_clip"),
    pivot="donor",
    figsize=(4, 4),
    vmin=0,
    **kwargs,
):
    for col in cols:
        sns.clustermap(col_corr(df, col, pivot=pivot), figsize=figsize, vmin=vmin, **kwargs).fig.suptitle(col)


def get_de_stats(df, metric_dict=frozendict({"": np.mean, "_var": np.var}), grp="gene"):
    stat_dfs = list()
    for suffix, metric in metric_dict.items():
        stat_dfs.append(df.groupby(grp).agg(metric).add_suffix(suffix))
    df = pd.concat(stat_dfs, axis=1)
    df["score"] = df["logfc"] * df["-logp"]
    df = df.sort_values(["score", "logfc"], ascending=False)
    return df
