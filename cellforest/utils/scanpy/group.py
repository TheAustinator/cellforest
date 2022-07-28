import os
from copy import copy
from pathlib import Path
from typing import Dict, Callable, Iterable, Union, Tuple, Optional, Literal, Any

import anndata
from anndata import AnnData
import numpy as np
import pandas as pd
from scanpy import read_h5ad

_AGG_FUNCS = {
    "sum": lambda x: x.sum(axis=0),
    "mean": lambda x: x.mean(axis=0),
    "var": lambda x: x.var(axis=0),
    "nonzero": lambda x: (x > 0).sum(axis=0),
}


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


def groupby_agg(ad, agg: Union[Literal["sum", "mean", "nonzero"], Callable], obs_cols: Iterable, layers: Optional[str] = None, obs_aggs: Optional[Dict[Any, Callable]] = None):
    agg = agg if callable(agg) else _AGG_FUNCS[agg]
    layers = [] if layers is None else layers
    ad_d = groupby_dict(ad, obs_cols)
    def _get_layer(_ad, lay): return _ad.layers[lay] if lay != "X" else ad.X
    ser_size = pd.Series({k: len(_ad) for k, _ad in ad_d.items()})
    ser_agg = {lay: pd.Series({k: agg(_get_layer(_ad, lay)) for k, _ad in ad_d.items()}) for lay in layers + ["X"]}
    obs_aggs_res = pd.DataFrame({col: {k: f(_ad.obs[col]) for k, _ad in ad_d.items()} for col, f in obs_aggs.items()})
    def _agg_layer(lay): return np.concatenate([np.expand_dims(x, 1) for x in ser_agg[lay].values.tolist()], axis=1).T
    obs = ser_agg["X"].index.to_frame().drop(columns=0).merge(obs_aggs_res, left_index=True, right_index=True, how="left")
    layers = {lay: _agg_layer(lay) for lay in layers}
    ad = AnnData(X=_agg_layer("X"), obs=obs, var=ad.var, layers=layers)
    ad.obs["size"] = ser_size
    return ad


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


def groupby_ctrl(
        ad: AnnData, obs_cols: Iterable, ctrl_col: str, ctrl_val: str, join: str = "outer", return_key: bool = True
):
    obs_cols = [obs_cols,] if isinstance(obs_cols, str) else list(obs_cols)
    ad_d = groupby_dict(ad, obs_cols)
    ind_i = obs_cols.index(ctrl_col)
    for k, _ad in ad_d.items():
        k_ctrl = list(copy(k))
        k_ctrl[ind_i] = ctrl_val
        k_ctrl = tuple(k_ctrl)
        if k != k_ctrl:
            _ad = anndata.concat([_ad, ad_d[k_ctrl]], join=join)
        if return_key:
            yield k, _ad
        else:
            yield _ad


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
