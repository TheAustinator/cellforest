import os
from copy import copy
from pathlib import Path
from typing import Dict, Callable, Iterable, Union, Tuple

import anndata
from anndata import AnnData
from scanpy import read_h5ad


def groupby(
    ad: AnnData, obs_cols: Iterable, return_key: bool = True
) -> Union[Tuple[str, AnnData], AnnData]:
    grp = ad.obs.groupby(obs_cols)
    for name, obs_sub in grp:
        selector = ad.obs_names.isin(obs_sub.index)
        _ad = ad[selector]
        _ad._sanitize()
        if return_key:
            yield name, _ad
        else:
            yield _ad


def groupby_dict(ad: AnnData, obs_cols: Iterable):
    ad_dict = dict()
    for name, ad_sub, in groupby(ad, obs_cols):
        ad_dict[name] = ad_sub
    return ad_dict


def groupby_dict_ctrl(
    ad: AnnData, obs_cols: Iterable, ctrl_col: str, ctrl_val: str, join: str = "outer",
):
    obs_cols = [obs_cols,] if isinstance(obs_cols, str) else list(obs_cols)
    ad_d = groupby_dict(ad, obs_cols)
    ind_i = obs_cols.index(ctrl_col)
    for k, _ad in ad_d.items():
        k_healthy = list(copy(k))
        k_healthy[ind_i] = ctrl_val
        k_healthy = tuple(k_healthy)
        if k != k_healthy:
            ad_d[k] = anndata.concat([_ad, ad_d[k_healthy]], join=join)
    return ad_d


def values_apply(
    d: Dict[str, AnnData], func: Callable, print_keys: bool = False, **kwargs
):
    output = dict()
    for k, v in d.items():
        if print_keys:
            print(k)
        output[k] = func(v, **kwargs)
    return output


def write(ad_d: Dict[str, AnnData], root_dir: str):
    os.makedirs(root_dir, exist_ok=True)
    for k, _ad in ad_d.items():
        _ad.write(Path(root_dir) / f"{k}.h5ad")


def read(root_dir: str):
    paths = Path(root_dir).glob("*.h5ad")
    return {p.stem: read_h5ad(p) for p in paths}
