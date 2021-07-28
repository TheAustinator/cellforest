from pathlib import Path
from typing import Optional

from anndata import AnnData
from frozendict import frozendict
import numpy as np
import pandas as pd

import seaborn as sns

from cellforest.api.r import ad_to_r
from cellforest.utils import r
from cellforest.utils.shell.shell_command import process_shell_command

_R_UTILS_DIR = Path(r.__file__).parent
_R_MAST_DISK = str(_R_UTILS_DIR / "_mast_disk.R")


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


def mast(ad: AnnData, formula: str, cores: int = -1, disk: bool = True, log_dir: str = "/tmp"):
    if disk:
        return _mast_disk(ad, formula, cores, log_dir)
    else:
        return _mast_mem(ad, formula)


def _mast_disk(ad: AnnData, formula: str, cores: int, log_dir: str = "/tmp"):
    ad_to_r(ad, "/tmp/ad_sce.rds", format="sce")
    cmd_str = f"Rscript {_R_MAST_DISK} {formula} {cores}"
    process_shell_command(cmd_str, log_dir, "mast")
    df = pd.read_csv("/tmp/df_zlm.csv", index_col=0)
    return df


def _mast_mem(
    ad: AnnData, formula: str,
):
    raise NotImplementedError("MAST in memory not yet supported.")
