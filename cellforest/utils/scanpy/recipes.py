from copy import deepcopy
from typing import List, Dict, Any, Optional

from dataforest.utils.structures.dict import map_dict_vals

from cellforest.utils.plot.diffexp import volcano
from cellforest.utils.scanpy.cell import markers, get_markers_df
from cellforest.utils.scanpy.group import groupby_dict


def map_list_dict_key(
    l: List[dict],
    k: Any,
    f: callable,
    new_key: Optional[str] = None,
    return_updated=True,
    key_kwargs: dict = None,
    **kwargs,
):
    """
    for a list of dicts, map a function to a given key. Returns a deepcopy of
    dict with new value at same key if `new_key=None`, otherwise, adds under
    `new_key` in original dict
    Args:
        l:
        k:
        f:
        new_key: new key in dict at which to place output in original dict. If
            not specified, deepcopy is returned with output at original key
        return_updated: returns nothing if False, but still updates in place if
            no `new_key`
        key_kwargs: mapping of {function_kwarg: dict_kwarg} to pass any kwargs
             to function from dict itself
        **kwargs:

    Returns:

    """
    # TODO: have a decorator for this somewhere -- replace
    key_kwargs = key_kwargs if key_kwargs else dict()
    if new_key is None and return_updated:
        l = deepcopy(l)
        new_key = k
    for d in l:
        _kwargs = {_k: d[_v] for _k, _v in key_kwargs.items()}
        _kwargs.update(kwargs)
        v = d[k]
        output = f(v, **_kwargs)
        if return_updated or new_key == k:
            d[new_key] = output
    if return_updated:
        return l


def flow_1(spec: List[Dict[str:Any]]):
    """

    Args:
        spec: [{"name": "diffexp_1", "ad": my_anndata, "key": "my_obs_col", "val": "my_obs_col_val"}, ...]

    Returns:

    """
    map_list_dict_key(spec, "ad", markers, return_updated=False, key_kwargs={"key": "key"})
    # group anndata by donor
    map_list_dict_key(spec, "ad", lambda ad: groupby_dict(ad, "donor"), new_key="ad_d")
    # TODO: broken
    # map_list_dict_key(spec, "ad_d", lambda d: map_dict_vals(d, lambda ad, _k: markers(ad)))
    map_list_dict_key(spec, "ad", get_markers_df, new_key="de")
    filter_de = lambda de, val: de[de["group"] == val]
    map_list_dict_key(spec, "de", filter_de, key_kwargs={"val": "val"})
    map_list_dict_key(spec, "de", volcano)
    map_list_dict_key(spec, "de", volcano, xlim=(-10, 10), pval_thresh=0.1, logfc_thresh=0.5)
