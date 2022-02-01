import gc
import logging
import os
from pathlib import Path
from typing import Optional, Dict, Union, Literal

from anndata import AnnData
from lazy_import import lazy_module, lazy_callable
anndata2ri = lazy_module("anndata2ri")
robjects = lazy_module("rpy2.robjects")
localconverter = lazy_callable("rpy2.robjects.conversion.localconverter")

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


def ad_to_r(ad: AnnData, filepath: AnyPath, format_: Literal["seurat", "sce"] = "sce", mem_gb=None):
    ad = prep(ad)
    with localconverter(anndata2ri.converter):
        ad_r = anndata2ri.py2rpy(ad)
    robjects.r.assign("ad_r", ad_r)
    if mem_gb:
        robjects.r(f"Sys.setenv('R_MAX_VSIZE'={mem_gb * 1000000000})")
    try:
        if format_ == "seurat":
            if "logcounts" not in ad.layers:
                raise ValueError('Seurat conversion requires "logcounts"')
            robjects.r("ad_r <- Seurat::as.Seurat(ad_r)")
            robjects.r("ad_r@assays$RNA <- ad_r@assays$originalexp")
            robjects.r("ad_r@assays$originalexp <- NULL")
        robjects.r(f"saveRDS(ad_r, file='{filepath}')")
        robjects.r("rm(ad_r)")
        robjects.r("gc()")
    except Exception as e:
        try:
            robjects.r("rm(ad_r)")
            robjects.r("gc()")
        except Exception:
            pass
        gc.collect()
        raise e
    gc.collect()


def clear_tmp(path: str = "/tmp", pattern: str = "./cf_*"):
    paths = [p for p in Path(path).glob(pattern)]
    for p in paths:
        os.remove(p)
