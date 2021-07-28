import gc

import logging
from pathlib import Path
from typing import Optional, Dict, Union, Literal

from anndata import AnnData
import anndata2ri
from rpy2 import robjects
from rpy2.robjects.conversion import localconverter

LAYER_RENAMES_DEFAULT = {"raw": "counts", "log": "logcounts"}

AnyPath = Union[str, Path]


def prep(ad: AnnData, layer_renames: Optional[Dict[str, str]] = None) -> AnnData:
    """
    For compatibility with anndata2ri re-instantiates and converts all object columns to cat
    """
    ad = AnnData(ad.X, ad.obs, ad.var, ad.uns, ad.obsm, ad.varm, ad.layers, ad.raw, obsp=ad.obsp, varp=ad.varp)
    layer_renames = layer_renames if layer_renames else LAYER_RENAMES_DEFAULT
    for layer, alias in layer_renames.items():
        if layer in ad.layers:
            ad.layers[alias] = ad.layers[layer]
        else:
            logging.warning(f"layer ({layer}) from layer renames missing -- skipping. `layer_renames={layer_renames}`")
    cols_obj = ad.obs.dtypes[ad.obs.dtypes == "object"].index.tolist()
    ad.obs[cols_obj] = ad.obs[cols_obj].astype("category")
    gc.collect()
    return ad


def ad_to_r(ad: AnnData, filepath: AnyPath, format: Literal["seurat", "sce"] = "sce"):
    ad = prep(ad)
    with localconverter(anndata2ri.converter):
        ad_r = anndata2ri.py2rpy(ad)
    robjects.r.assign("ad_r", ad_r)
    try:
        if format == "seurat":
            robjects.r("ad_r <- Seurat::as.Seurat(ad_r)")
        robjects.r(f"saveRDS(ad_r, file='{filepath}')")
        robjects.r("gc()")
    except Exception as e:
        try:
            robjects.r("gc()")
        except:
            pass
        gc.collect()
        raise e
    gc.collect()
