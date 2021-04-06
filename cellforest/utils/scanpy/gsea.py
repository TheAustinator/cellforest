from typing import Union, AnyStr, List, Dict, Optional, Callable

from anndata import AnnData
from dataforest.utils.analysis import set as set_metrics
from dataforest.utils.analysis.pairwise import pairwise_metric
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from rpy2 import robjects

from cellforest.utils import parse_gene_set_gmt
from cellforest.utils.scanpy.gene import rank_markers
from cellforest.utils.scanpy.manifold import pc_load_df

R_GSEA_STR = """
r_gsea <- function(genes, ranks, gmt_path) {
    set.seed(54321)
    library(data.table)
    library(fgsea)
    library(gage)
    library(ggplot2)
    
    names(ranks) <- genes
    ranks <- unlist(ranks)
    
    gene_set = fgsea::gmtPathways(gmt_path)
    df_fgsea <- fgsea(pathways = gene_set,
                  stats    = ranks,
                  minSize  = 15,
                  maxSize  = 500)
    df_fgsea <- as.data.frame(df_fgsea)
    df_fgsea$leadingEdge <- sapply(df_fgsea$leadingEdge, toString)
    
    df_gage <- gage::gage(ranks, gsets=gene_set, same.dir=TRUE, set.size =c(15,600))
    
    write.csv(df_fgsea, "fgsea_test.csv")
    write.csv(df_gage, "gage_test.csv")
}
"""


def add_ranks(df, rank_max_fc=0.5):
    def _update_saturated_ranks(df_sub, extrema):
        scale_factor = 1 + rank_max_fc * df_sub["logfc"].abs() / (df_sub["logfc"].max() - df_sub["logfc"].min())
        df_sub["rank"] = extrema
        df_sub["rank"] *= scale_factor
        return df_sub

    df["rank"] = rank_markers(df)
    sorted_ranks = sorted(df["rank"].unique())
    rank_min = sorted_ranks[1]
    rank_max = sorted_ranks[-2]
    max_selector = df["rank"] == np.inf
    min_selector = df["rank"] == -np.inf
    df.loc[max_selector] = _update_saturated_ranks(df.loc[max_selector], rank_max)
    df.loc[min_selector] = _update_saturated_ranks(df.loc[min_selector], rank_min)
    df.sort_values(["rank", "logfc"], inplace=True)
    genes = df["gene"].tolist()
    ranks = df["rank"].tolist()
    return genes, ranks


def gsea(df, gmt_path, ranked=False, rank_max_fc=0.5):
    """

    Args:
        df:
        gmt_path:
        ranked: whether genes are already ranked. If so, expected to be ordered with
            a column, "ranked", which contains some sort of cardinal or ordinal score.
            Otherwise, expecting "logfc" and "pval_adj" columns for ranking.
        rank_max_fc: modulates dynamic range of ranking metric for values where logp=inf

    Returns:

    """
    r_gsea = robjects.r(R_GSEA_STR)
    if ranked:
        genes = df.index.tolist()
        ranks = df["rank"].tolist()
    else:
        genes, ranks = add_ranks(df, rank_max_fc)
    r_gsea(genes, ranks, str(gmt_path))
    df_fgsea = pd.read_csv("fgsea_test.csv", index_col=0)
    df_fgsea.set_index("pathway", inplace=True)
    df_gage = pd.read_csv("gage_test.csv", index_col=0)
    df_gsea = df_fgsea.merge(df_gage, left_index=True, right_index=True)
    return df_gsea

    # numerics = ['int16', 'int32', 'int64', 'float16', 'float32', 'float64']
    # df_gsea_num = df_gsea.select_dtypes(include=numerics)
    # corr = df_gsea_num.corr()
    # num_fgsea = df_fgsea.select_dtypes(include=numerics).columns
    # num_gage = df_gage.select_dtypes(include=numerics).columns
    # corr.loc[num_fgsea, num_fgsea] = -1
    # corr.loc[num_gage, num_gage] = -1
    # clustermap(corr)
    # plt.scatter(df_gsea["stats.stat.mean"], df_gsea["NES"], s=1, alpha=0.1)
    # plt.scatter(df_gsea["greater.q.val"], df_gsea["padj"], s=1, alpha=0.1)


def pc_gsea(
    ad: AnnData, gmt_path: AnyStr, pval: float = 0.01, return_df=True, parallel=False,
) -> Union[Dict[str, List[str]], pd.DataFrame]:
    pc_df = pc_load_df(ad)

    def _get_top_gene_sets(series: pd.Series):
        _df_gsea = gsea(pd.DataFrame({"rank": series}), gmt_path, ranked=True)
        top_gs = _df_gsea[_df_gsea["padj"] < pval].sort_values("NES", ascending=False).index.tolist()
        return top_gs

    if parallel:
        raise NotImplementedError("Parallel doesn't work with rpy2 b/c stateful")
        # top_gs_list = Parallel(n_jobs=-1)(delayed(_get_top_gene_sets)(pc_df[col]) for col in pc_df)
    else:
        top_gs_list = [_get_top_gene_sets(pc_df[col]) for col in pc_df]
    pc_gs = dict(zip(pc_df.columns.tolist(), top_gs_list))
    if return_df:
        max_len = max(map(len, pc_gs.values()))
        pc_gs = {k: v + (max_len - len(v)) * [""] for k, v in pc_gs.items()}
        pc_gs = pd.DataFrame(pc_gs)
    return pc_gs


def gene_set_pair_metric(
    gmt: Union[dict, AnyStr], gene_sets: Optional[List[str]] = None, metric: Union[AnyStr, Callable] = "iou"
):
    if not isinstance(gmt, dict):
        gmt = parse_gene_set_gmt(gmt)
    if gene_sets:
        gmt = {k: v for k, v in gmt.items() if k in gene_sets}
    if isinstance(metric, str):
        metric = getattr(set_metrics, metric)
    df = pairwise_metric(gmt, metric)
    return df
