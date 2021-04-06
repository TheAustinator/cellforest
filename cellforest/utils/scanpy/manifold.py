from typing import Callable, Optional

from anndata import AnnData
from dataforest.utils.structures.dict import map_dict_vals
import pandas as pd

from cellforest.utils.scanpy.cell import reduce
from cellforest.utils.scanpy.group import groupby_dict


def pc_load_df(ad: AnnData, n_pcs=None, drop_zero=True):
    X = ad.varm["PCs"]
    n_pcs = n_pcs if n_pcs else X.shape[1]
    X = X[:, :n_pcs]
    cols = [f"PC_{i}" for i in range(n_pcs)]
    df = pd.DataFrame(X, index=ad.var_names, columns=cols)
    if drop_zero:
        df = df.loc[df.sum(axis=1) != 0]
    return df


def pc_load_df_grp(
    ad: AnnData, grp_key: str, preprocess_func: Optional[Callable] = reduce, n_pcs: Optional[int] = None
):
    """

    Args:
        ad:
        grp_key: `ad.obs` key for grouping
        preprocess_func: function to preprocess sub-grouped anndatas
        n_pcs:

    Returns:

    """
    ad_d = groupby_dict(ad, grp_key)
    if preprocess_func:
        map_dict_vals(ad_d, preprocess_func)
    pc_df_dict = map_dict_vals(
        ad_d, lambda key, _ad: pc_load_df(_ad, n_pcs, drop_zero=False).add_prefix(f"{key}_"), pass_keys=True
    )
    return pd.concat(list(pc_df_dict.values()), axis=1)
