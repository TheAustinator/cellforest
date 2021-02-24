from anndata import AnnData
from typing import Dict, Callable, Iterable


def groupby(ad: AnnData, obs_cols: Iterable):
    grp = ad.obs.groupby(obs_cols)
    for name, obs_sub in grp:
        selector = ad.obs_names.isin(obs_sub.index)
        yield name, ad[selector]


def groupby_dict(ad: AnnData, obs_cols: Iterable):
    ad_dict = dict()
    for name, ad_sub, in groupby(ad, obs_cols):
        ad_dict[name] = ad_sub
    return ad_dict


def values_apply(d: Dict[str, AnnData], func: Callable, **kwargs):
    output = dict()
    for k, v in d.items():
        output[k] = func(v, **kwargs)
    return output
