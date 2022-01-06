from functools import wraps
from typing import Optional, Dict, Tuple, Union, Set, List

from anndata import AnnData
from dataforest.utils.decorators import sub_kwargs
import numpy as np
import pandas as pd
import scanpy as sc

from cellforest.utils.scanpy.generic import _generic_preprocess

MARKER_THRESH_DEFAULTS = {
    "gt": {"logfc": 0.25},
    "lt": {"pval_adj": 0.05},
}


def _check_mat_type(X) -> str:
    sums = X.sum(axis=1)
    sums = pd.Series(np.array(sums).flatten())
    is_raw = sums.map(float.is_integer).mean() > 0.9
    if is_raw:
        return "raw"
    is_norm = sums.var() / sums.mean() < 1
    if is_norm:
        return "norm"
    if (sums > 0).mean() < 0.9:
        return "scale"
    else:
        return "log"


def preprocess(
    ad: AnnData,
    min_cells: int = None,
    min_genes: int = None,
    max_genes: int = None,
    max_pct_mito: int = None,
):
    x_type = _check_mat_type(X)

    if "raw" not in ad.layers:
        ad.layers["raw"] = ad.X.copy()
    ad.raw = ad.copy()
    ad = _generic_preprocess(ad, min_cells, min_genes, max_genes, max_pct_mito)
    ad.layers["norm"] = ad.X.copy()
    sc.pp.log1p(ad)
    ad.layers["log"] = ad.X
    return ad


@sub_kwargs(
    "hvg_kwargs", "pca_kwargs", "neighbors_kwargs", "umap_kwargs", "tsne_kwargs"
)
def reduce(
    ad: AnnData,
    n_comps: int = 15,
    hvg_kwargs=None,
    pca_kwargs=None,
    neighbors_kwargs=None,
    umap_kwargs=None,
    tsne_kwargs=None,
    umap=True,
    tsne=False,
    scale: bool = False,
):
    sc.pp.highly_variable_genes(ad, batch_key="sample", **hvg_kwargs)
    if scale:
        sc.pp.scale(ad)
    sc.pp.pca(ad, n_comps=n_comps, **pca_kwargs)
    sc.pp.neighbors(ad, **neighbors_kwargs)
    if umap:
        sc.tl.umap(ad, **umap_kwargs)
    if tsne:
        sc.tl.tsne(ad, **tsne_kwargs)
    return ad


def cluster(
    ad: AnnData,
    neighbors_key: str = None,
    embeddings_key: str = "X_umap",
    resolution: float = 0.1,
    do_markers=True,
    key_added=None,
    **markers_kwargs,
):
    if neighbors_key and not key_added:
        key_added = f"leiden_{neighbors_key}"
    elif not key_added:
        key_added = "leiden"
    sc.tl.leiden(
        ad, neighbors_key=neighbors_key, resolution=resolution, key_added=key_added
    )
    if do_markers:
        markers(ad, key=key_added, embeddings_key=embeddings_key, **markers_kwargs)


def markers(
    ad: AnnData,
    key: str = "leiden",
    embeddings_key: str = "X_umap",
    layer: str = "log",
    n_genes: int = 30,
    filter_var_prefix=(),
    obs_facet=None,
    fontsize: int = 8,
    ax=None,
    plot_clusters=True,
    plot_markers=True,
):
    if filter_var_prefix:
        ad = ad[:, ~ad.var_names.str.startswith(filter_var_prefix)]
    if obs_facet:
        for val, _obs in ad.obs.groupby(obs_facet):
            print(val)
            _ad = ad[ad.obs_names.isin(_obs.index)]
            if plot_clusters:
                sc.pl.embedding(_ad, embeddings_key, color=key)
    else:
        if plot_clusters:
            sc.pl.embedding(ad, embeddings_key, color=key)
    sc.tl.rank_genes_groups(ad, groupby=key, layer=layer, use_raw=False)
    if plot_markers:
        sc.pl.rank_genes_groups(ad, n_genes=n_genes, fontsize=fontsize, ax=ax)
    get_markers_df(ad)


def marker_dict(
    ad,
    gt_thresh: Optional[Dict[str, float]] = None,
    lt_thresh: Optional[Dict[str, float]] = None,
    sort_by: str = "volcano",
    sort_ascending: bool = False,
    max_n: Optional[int] = None,
    prefix_filter: Optional[Union[List[str], str]] = None,
    return_stats: bool = True,
) -> Union[Tuple[Dict[str, List[str]], pd.DataFrame], Dict[str, List[str]]]:
    gt_thresh = gt_thresh if gt_thresh is not None else MARKER_THRESH_DEFAULTS["gt"]
    lt_thresh = lt_thresh if lt_thresh is not None else MARKER_THRESH_DEFAULTS["lt"]
    mark = ad.uns["markers"]
    mark.sort_values(sort_by, ascending=sort_ascending, inplace=True)
    genes_dict = dict()
    stats_dict = dict()
    for name, markers_grp in mark.groupby("group"):
        for col, thresh in gt_thresh.items():
            markers_grp = markers_grp[markers_grp[col] > thresh]
        for col, thresh in lt_thresh.items():
            markers_grp = markers_grp[markers_grp[col] < thresh]
        markers_grp["ribo"] = markers_grp["gene"].str.startswith(("RPL", "RPS"))
        markers_grp["hsp"] = markers_grp["gene"].str.startswith("HSP")
        markers_grp["mt"] = markers_grp["gene"].str.startswith("MT-")
        stats = {
            "max_pval": markers_grp["pval_adj"].max(),
            "min_logfc": markers_grp["logfc"].min(),
            "n_markers": len(markers_grp),
            "frac_ribo": markers_grp["ribo"].sum() / len(markers_grp),
            "frac_hsp": markers_grp["hsp"].sum() / len(markers_grp),
            "frac_mt": markers_grp["mt"].sum() / len(markers_grp),
        }
        genes_dict[name] = markers_grp["gene"].tolist()
        stats_dict[name] = stats
    prefix_filter = (
        [prefix_filter,] if isinstance(prefix_filter, str) else prefix_filter
    )
    for k, gene_list in genes_dict.items():
        if prefix_filter is not None:
            gene_list = [
                x
                for x in gene_list
                if not any([x.startswith(pre) for pre in prefix_filter])
            ]
        if max_n:
            gene_list = gene_list[: min(max_n, len(gene_list))]
        genes_dict[k] = gene_list
    df_stats = pd.DataFrame(stats_dict).T
    ret = (genes_dict, df_stats) if return_stats else genes_dict
    return ret


def process(ad: AnnData):
    ad = preprocess(ad)
    ad = reduce(ad)
    return ad


def get_markers_df(
    ad: AnnData, mean_expr_clip: float = 0.1, group=None, uns_key="markers"
):
    """

    Args:
        ad:
        mean_expr_clip:
        group:

    Returns:

    """

    def get_col(col):
        return np.array(ad.uns["rank_genes_groups"][col].tolist())

    cols = ["names", "pvals_adj", "logfoldchanges"]
    col_dict = {col: get_col(col).flatten() for col in cols}
    labels = ad.uns["rank_genes_groups"]["names"].dtype.names
    size = len(get_col("names"))
    groups = np.tile(labels, (size, 1)).flatten()
    col_dict["group"] = groups
    df = pd.DataFrame(col_dict)
    df.rename(
        columns={"names": "gene", "pvals_adj": "pval_adj", "logfoldchanges": "logfc"},
        inplace=True,
    )
    group_var = ad.uns["rank_genes_groups"]["params"]["groupby"]

    def _add_group(_group, _dff):
        selector = ad.obs[group_var] == _group
        if selector.sum() == 0:
            raise ValueError(
                f"specified `group`: {group}   not in column :{group_var}  . Choose from {ad.obs[group_var].unique()}"
            )
        mean_expr = ad[ad.obs[group_var] == _group, _dff["gene"]].X.mean(axis=0)
        frac_expr = (ad[ad.obs[group_var] == _group, _dff["gene"]].X > 0).mean(axis=0)
        df.loc[_dff.index, "mean_expr"] = mean_expr.T
        df.loc[_dff.index, "frac_expr"] = frac_expr.T

    if not group:
        for group, dff in df.groupby("group"):
            _add_group(group, dff)
    else:
        _add_group(group, df)
    df["mean_expr_clip"] = df["mean_expr"].apply(
        lambda x: np.clip(x, 0, mean_expr_clip)
    )
    df["-logp"] = -df["pval_adj"].apply(np.log10)
    ceil = df["-logp"][df["-logp"].abs() != np.inf].max()
    df["-logp"] = df["-logp"].replace({np.inf: ceil})
    df["volcano"] = df["-logp"] * df["logfc"]
    df["volcano_abs"] = df["volcano"].abs()
    if uns_key:
        ad.uns[uns_key] = df
    return df


def annotate(
    ad: AnnData,
    col: str,
    basis_col: str,
    val: Union[str, dict],
    basis_vals: Optional[Union[str, List[str]]] = None,
):
    def _single_annotate(_ad: AnnData, _col: str, _basis_col: str, _val, _basis_vals):
        _basis_vals = (
            [_basis_vals,]
            if isinstance(_basis_vals, (int, float, str))
            else _basis_vals
        )
        ad.obs.loc[ad.obs[_basis_col].isin(_basis_vals), _col] = _val

    if col in ad.obs.columns:
        ad.obs[col] = ad.obs[col].astype(str)
    if isinstance(val, dict):
        for k, v in val.items():
            _single_annotate(ad, col, basis_col, v, k)
    else:
        _single_annotate(ad, col, basis_col, val, basis_vals)
    ad.obs[col] = ad.obs[col].astype("category")
