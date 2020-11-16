from pathlib import Path
from typing import Optional, Dict, Set, Union, Literal, List

import pandas as pd
from dataforest.utils import listify
from seaborn import barplot

from cellforest.utils import parse_gene_set_gmt


def gene_set_bar(
    genes: Set[str],
    gene_set_dict: Optional[Dict[str, Set[str]]] = None,
    gene_set_gmt_path: Optional[Union[str, Path]] = None,
    metric: Union[List[str], Literal["genes", "gene_set", "iou", "all"]] = "all",
    min_frac: float = 0.05,
    xtick_label_rot: int = 90,
    **kwargs,
):
    """
    Args:
        genes:
        gene_set_dict:
        gene_set_gmt_path:
        metric: str or list of str
            "genes" -- show fraction of your `genes` contained in each gene set;
            "gene_set" -- show fraction of each gene set contained in your `genes`;
            "iou" -- intersection over union of `genes` and gene set;
            "all" -- all three metrics as different colors;
        min_frac: min fraction to show in plot
    """

    def _io1(gs_1: Set[str], gs_2: Set[str]):
        return len(gs_1.intersection(gs_2)) / len(gs_1)

    def _iou(gs_1: Set[str], gs_2: Set[str]):
        return len(gs_1.intersection(gs_2)) / len(gs_1.union(gs_2))

    metric_func_lookup = {
        "genes": lambda _genes, _gene_set: _io1(genes, gene_set),
        "gene_set": lambda _genes, _gene_set: _io1(gene_set, genes),
        "iou": lambda _genes, _gene_set: _iou(genes, gene_set),
    }

    genes = set(genes)
    if not gene_set_dict or gene_set_gmt_path:
        raise ValueError("Must provide `gene_set_dict` or `gene_set_gmt_path`")
    elif gene_set_gmt_path:
        gene_set_dict = parse_gene_set_gmt(gene_set_gmt_path)
    metrics = listify(metric) if metric != "all" else ["genes", "gene_set", "iou"]
    rows = list()
    for name, gene_set in gene_set_dict.items():
        for m in metrics:
            func = metric_func_lookup[m]
            gene_set_metric = {
                "gene_set": name,
                "metric": m,
                "fraction": func(genes, gene_set),
            }
            rows.append(gene_set_metric)
    df = pd.DataFrame(rows)
    df.sort_values("fraction", inplace=True, ascending=False)
    df = df[df["fraction"] > min_frac]
    ax = barplot(data=df, x="gene_set", y="fraction", hue="metric", **kwargs)
    ax.xaxis.set_tick_params(rotation=xtick_label_rot)
    ax.axhline(min_frac)
    return ax
