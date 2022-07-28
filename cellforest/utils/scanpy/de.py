import gc
import logging
import os
from itertools import product
from pathlib import Path
from typing import Optional, Union, Iterable, Dict
import uuid

from anndata import AnnData
from frozendict import frozendict
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

from cellforest.api.r import ad_to_r
from cellforest.api.tl import undersample
from cellforest.utils import r
from cellforest.utils.scanpy.cell import get_markers_df
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
    ad: AnnData, formula: str, cores: int = 0, disk: bool = True, log_dir: str = "/tmp", keep_temp: bool = False,
) -> pd.DataFrame:
    ad = ad.copy()
    ad.obs = _str_dtypes(ad.obs)
    if disk:
        return _mast_disk(ad, formula, cores, log_dir, keep_temp)
    else:
        return _mast_mem(ad, formula)


def pairwise(ad: AnnData, obs_key: str = "sample"):
    de = pd.DataFrame()
    df_meta = pd.DataFrame()
    pairs = list(product(ad.obs[obs_key].unique(), ad.obs[obs_key].unique()))
    for pair in pairs:
        if pair[0] == pair[1]:
            continue
        ad.obs["grp"] = ""
        obs = ad.obs[ad.obs[obs_key].isin(pair)]
        counts = obs[obs_key].value_counts()
        n_cells = counts[counts > 0].min()
        df_meta[pair] = pd.Series({"n_cells": n_cells})
        grp_ids = {
            k: list(np.random.choice(_obs.index.tolist(), n_cells, replace=False))
            if len(_obs) > 0 else []
            for k, _obs in obs.groupby(obs_key)
        }
        for grp, ids in grp_ids.items():
            ad.obs.loc[ad.obs_names.isin(ids), "grp"] = grp
        ad.obs["grp"] = ad.obs["grp"].astype("category")
        sc.tl.rank_genes_groups(ad, "grp", groups=pair, reference=pair[1])
        _de = get_markers_df(ad).set_index("gene")
        for grp, dff in _de.groupby("group"):
            de[pair[0], pair[1], grp, "-logp"] = dff["-logp"]
            de[pair[0], pair[1], grp, "logfc"] = dff["logfc"]
            de[pair[0], pair[1], grp, "volcano"] = dff["volcano"]
    de = de.T
    df_meta = df_meta.T
    de.index = pd.MultiIndex.from_tuples(de.index)
    return de, df_meta


def grouped(ad, obs_cmp: str, obs_groupby: str = "sample", reference: Optional = None, test: Optional = None, undersample_bootstrap: bool = False, min_cells_skip: int = 10, as_anndata: bool = True) -> pd.DataFrame:
    de = pd.DataFrame()
    if undersample_bootstrap not in [False, 1]:
        raise NotImplementedError("bootstrapping not yet implemented, so must be False or 1")
    for grp in ad.obs[obs_groupby].unique():
        print(grp)
        ad.obs["_grp"] = np.nan
        grp_mask = ad.obs[obs_groupby] == grp
        ad.obs.loc[grp_mask, "_grp"] = ad.obs.loc[grp_mask, obs_cmp]
        ad.obs["_grp"] = ad.obs["_grp"].astype("category")
        def _notna(x): np.isnan(x) if isinstance(x, float) else False
        groups = sorted([x for x in ad.obs.loc[grp_mask, "_grp"].unique() if not _notna(x)])
        test = test if test is not None else set(groups).difference({reference}).pop()
        if len(groups) != 2:
            logging.info(f"len groups for {grp}: {groups} is not two, skiping.")
            continue
        min_ = ad.obs["_grp"].value_counts().min()
        if min_ < min_cells_skip:
            logging.info(f"minimum class is {min_} for one group in {grp}. skipping")
            continue
        if undersample_bootstrap:
            max_class = ad.obs["_grp"].value_counts().min()
            inds = list()
            for class_ in ad.obs["_grp"].unique():
                _inds = ad.obs.reset_index()[ad.obs.reset_index()["_grp"] == class_].index
                inds += np.random.choice(_inds, min(len(_inds), max_class), replace=False).tolist()
            mask = set(ad.obs_names).difference(ad.obs_names[inds])
            ad.obs.loc[mask, "_grp"] = np.nan

        reference = reference if reference is not None else groups[1]
        sc.tl.rank_genes_groups(ad, groupby="_grp", reference=reference)
        _de = get_markers_df(ad, group=test, skip_mean_expr=True)
        _de = _de[_de["group"] == str(test)].set_index("gene")
        # import ipdb; ipdb.set_trace()
        de[grp, "-logp"] = _de["-logp"]
        de[grp, "logfc"] = _de["logfc"]
        de[grp, "volcano"] = _de["volcano"]

    de = de.T
    de.index = pd.MultiIndex.from_tuples(de.index)
    if as_anndata:
        layers = {k: dff.drop(columns="level_1") for k, dff in de.reset_index(1).groupby("level_1")}
        X = layers["volcano"]
        de = AnnData(X, layers=layers, var=X.columns.to_frame().drop(columns="gene"), obs=X.index.to_frame().drop(columns=0))
    return de


def pairwise_grouped(ad: AnnData, obs_key: str = "sample", obs_groupby: str = "cell_type"):
    raise NotImplementedError()


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
