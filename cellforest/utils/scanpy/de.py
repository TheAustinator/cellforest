import gc
import logging
import os
from pathlib import Path
from typing import Optional, Union, Iterable, Dict
import uuid

from anndata import AnnData
from frozendict import frozendict
import numpy as np
import pandas as pd

import seaborn as sns

from cellforest.api.r import ad_to_r
from cellforest.api.tl import undersample
from cellforest.utils import r
from cellforest.utils.scanpy.group import groupby_dict, values_apply
from cellforest.utils.shell.shell_command import process_shell_command

_R_UTILS_DIR = Path(r.__file__).parent
_R_MAST_DISK = str(_R_UTILS_DIR / "_mast_disk.R")


def _str_dtypes(df):
    return df.loc[:, (df.dtypes == "object") | (df.dtypes == "category")].astype(str)


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
        sns.clustermap(
            col_corr(df, col, pivot=pivot), figsize=figsize, vmin=vmin, **kwargs
        ).fig.suptitle(col)


def get_de_stats(df, metric_dict=frozendict({"": np.mean, "_var": np.var}), grp="gene"):
    stat_dfs = list()
    for suffix, metric in metric_dict.items():
        stat_dfs.append(df.groupby(grp).agg(metric).add_suffix(suffix))
    df = pd.concat(stat_dfs, axis=1)
    df["score"] = df["logfc"] * df["-logp"]
    df = df.sort_values(["score", "logfc"], ascending=False)
    return df


def mast(
    ad: AnnData, formula: str, cores: int = 0, disk: bool = True, log_dir: str = "/tmp"
) -> pd.DataFrame:
    ad = ad.copy()
    ad.obs = _str_dtypes(ad.obs)
    if disk:
        return _mast_disk(ad, formula, cores, log_dir)
    else:
        return _mast_mem(ad, formula)


def mast_grouped(
    ad_d: Union[AnnData, Dict[str, AnnData]],
    formula: str,
    grp: Optional[Union[Iterable[str], str]],
    undersample_col: Optional[str] = None,
    undersample_class_max: Optional[int] = None,
    cores: int = 0,
    disk: bool = True,
    log_dir: str = "/tmp",
    skip_error: bool = True,
    checkpoint_save_path: Optional[str] = "/tmp/mast_checkpoint.csv",
    resume: bool = False,
):
    if isinstance(ad_d, AnnData):
        ad_d = groupby_dict(ad_d, grp)

    if undersample_col:
        sampler = lambda _ad: undersample(_ad, undersample_col, undersample_class_max)
        ad_d = values_apply(ad_d, sampler)
    if resume and Path(checkpoint_save_path).exists():
        de = pd.read_csv(checkpoint_save_path)
        done = set(list(map(tuple, de[list(grp)].to_records(index=False))))
    else:
        de = pd.DataFrame()
        done = set()
    for k, _ad in ad_d.items():
        print(k)
        if k in done:
            print(f"Exists in checkpoint and `resume=True`. Skipping {k}")
            if skip_error:
                continue
            else:
                raise e
        if undersample_col:
            print(_ad.obs[undersample_col].value_counts())
        try:
            _de = mast(_ad, formula, cores, disk, log_dir)
        except Exception as e:
            print(f"Error on {k}: {e}")
            if skip_error:
                continue
            else:
                raise e
        if not isinstance(_de, pd.DataFrame):
            print(f"Error on {k}: {_de}")
        for (g, _k) in zip(grp, k):
            _de[g] = _k
        de = pd.concat([de, _de])
        if checkpoint_save_path:
            de.to_csv(checkpoint_save_path)
    return de


def _mast_disk(
    ad: AnnData, formula: str, cores: int, log_dir: str = "/tmp", keep_tmp: bool = False
):
    run_id = str(uuid.uuid4())[:4]
    sce_path = f"/tmp/cf_ad_sce_{run_id}.rds"
    formula = '"' + formula.replace(" ", "") + '"'
    ad_to_r(ad, sce_path, format_="sce")
    cmd_str = f"Rscript {_R_MAST_DISK} {run_id} {formula} {cores}"
    logging.info(f"Running MAST with run ID: {run_id}")
    logging.info(f"com")
    try:
        process_shell_command(cmd_str, log_dir, f"cf_mast_{run_id}")
    except Exception as e:
        _rm_tmp(sce_path, keep_tmp)
        raise e
    else:
        _rm_tmp(sce_path, keep_tmp)
    df = pd.read_csv(f"/tmp/cf_df_zlm_{run_id}.csv", index_col=0)
    gc.collect()
    return df


def _mast_mem(
    ad: AnnData, formula: str,
):
    raise NotImplementedError(
        "MAST in memory not yet supported. Parallelism didn't work via rpy2."
    )


def _rm_tmp(path: str, keep: bool):
    if not keep:
        os.remove(path)
