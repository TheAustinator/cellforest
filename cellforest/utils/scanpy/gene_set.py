import logging
from pathlib import Path
from typing import List, Union, Dict, AnyStr, Set, Literal, Optional

from anndata import AnnData
from dataforest.utils.analysis.set import iom, iou
import pandas as pd

from cellforest.utils import parse_gene_set_gmt
from cellforest.utils.scanpy.plot import embedding


def query_gene_sets(genes: Union[str, List[str]], gmt: Union[AnyStr, Dict[str, str]], iom_min: float = 1.0):
    if not isinstance(gmt, dict):
        gmt = parse_gene_set_gmt(gmt)
    gs_list = [
        {"gene_set": k, "iom": iom(genes, gs), "iou": iou(genes, gs), "len": len(gs)}
        for k, gs in gmt.items()
        if iom(genes, gs) >= iom_min
    ]
    return pd.DataFrame(gs_list)


def add_sigs(
    ad,
    gmt_dict: Optional[Dict[str, Set[str]]] = None,
    gene_sets: Optional[List[str]] = None,
    vmax_frac=0.95,
    plot=True,
    **kwargs,
):
    if not gmt_dict:
        raise NotImplementedError("No default gmt set up yet")
    if isinstance(gmt_dict, (str, Path)):
        gmt_dict = parse_gene_set_gmt(gmt_dict)
    if not gene_sets:
        gene_sets = list(gmt_dict.keys())
    gene_sets = [gene_sets,] if isinstance(gene_sets, str) else gene_sets
    gene_sets = list(gene_sets)
    _gene_sets = list()
    for set_name in gene_sets:
        genes = gmt_dict[set_name]
        if set(genes).intersection(ad.var_names):
            _gene_sets.append(set_name)
            add_gene_set(ad, genes, set_name, plot=False, **kwargs)
    if plot:
        embedding(ad, color=_gene_sets, vmax_frac=vmax_frac)


def add_gene_set(ad, genes, name, binary=False, gene_norm=False, plot=True, verbose=True):
    missing = [g for g in genes if g not in ad.var_names]
    union = [g for g in genes if g in ad.var_names]
    if missing:
        logging.warning(f"{len(union)}/{len(genes)} found from {name}. Remainder skipped.")
        if verbose:
            logging.warning(f"missing genes: {missing}")
        if not union:
            raise ValueError(f"No genes from {name} present")
    X = ad.raw[:, list(union)].X
    if binary:
        X = X > 0
    elif gene_norm:
        X = ad.X[:, ad.X.toarray().sum(axis=0) > 0]
        X = X / X.sum(axis=0)
    ad.obs[name] = X.sum(axis=1)
    if plot:
        embedding(ad, color=name, vmax_frac=0.9)


def gene_sets_union(gs_names: List[str], gmt: Union[str, Dict[str, Set[str]]]) -> Set[str]:
    if not isinstance(gmt, dict):
        gmt = parse_gene_set_gmt(gmt)
    return set.union(*[gmt[gs] for gs in gs_names])


def cluster_gene_rankings(ad: AnnData, mode: Literal["pval", "logfc", "volcano"] = "volcano") -> pd.DataFrame:
    mode_col_lookup = {"pval": "pval_adj"}
    mode_col = mode_col_lookup.get(mode, mode)
    df = pd.DataFrame({name: dff for name, dff in ad.uns["markers"].set_index("gene").groupby("group")[mode_col]})
    return df
