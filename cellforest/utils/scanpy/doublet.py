import logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc


SOLO_COLS = ["solo_pred", "solo_dub", "solo_score"]
DUB_FIND_COLS = ["dub_find_pred"]
SCRUBLET_COLS = ["scrublet_pred"]
ALL_COLS = SOLO_COLS + DUB_FIND_COLS + SCRUBLET_COLS
DOUBLET_LABELS = {
    "solo_pred": True,
    "solo_dub": True,
    "dub_find_pred": "Doublet",
}


def add_solo(ad, root_path):
    for label in ad.obs[["sample", "cell_types"]].agg("-".join, axis=1).unique():
        path = Path(f"{root_path}/solo/{label}/outputs.tsv")
        try:
            df = pd.read_csv(path, sep="\t", index_col=0)
        except FileNotFoundError:
            logging.warning(f"solo output not found for {path}. Skipping")
            ad.obs[SOLO_COLS] = np.nan
        else:
            df = df[df.index.isin(ad.obs.index)]
            ad.obs.loc[df.index, SOLO_COLS] = df
    return ad


def add_dub_find(ad, root_path):
    for label in ad.obs[["sample", "cell_types"]].agg("-".join, axis=1).unique():
        path = Path(f"{root_path}/doublet_finder/{label}/output.csv")
        try:
            df = pd.read_csv(path, sep=" ", index_col=0, header=None)[1]
            df = df.replace({"Singlet": False, "Doublet": True})
        except FileNotFoundError:
            logging.warning(f"DoubletFinder output not found for {path}. Skipping")
            ad.obs["dub_find_pred"] = np.nan
        else:
            df = df[df.index.isin(ad.obs.index)]
            ad.obs.loc[df.index, "dub_find_pred"] = df
    return ad


def _add_dub_methods(ad, root_path):
    if root_path is None:
        raise ValueError(f"`root_path` must be specified to locate doublet results")
    if "solo_pred" not in ad.obs.columns:
        ad = add_solo(ad, root_path)
    if "dub_find_pred" not in ad.obs.columns:
        ad = add_dub_find(ad, root_path)
    return ad


def lane_doublet_rates(ad, root_path=None, est_per_1k=0.008):
    ad = _add_dub_methods(ad, root_path)
    doublet_rates = []
    for label in ad.obs[["sample", "cell_types"]].agg("-".join, axis=1).unique():
        sample, cell_types = label.split("-")
        _ad = ad[(ad.obs["sample"] == sample) & (ad.obs["cell_types"] == cell_types)]
        # _ad.obs["dub_find_pred"] = _ad.obs["dub_find_pred"].replace({"Singlet": False, "Doublet": True})
        est = est_per_1k * len(_ad) / 1000
        solo_pred = (_ad.obs["solo_pred"] == "True").sum() / len(_ad)
        solo_dub = (_ad.obs["solo_dub"] == "True").sum() / len(_ad)
        dub_find_pred = (_ad.obs["dub_find_pred"] == "True").sum() / len(_ad)
        # scrublet_pred = _ad.obs["scrublet_pred"].sum() / len(_ad)
        row = {
            "label": label,
            "est": est,
            "solo_pred": solo_pred,
            "solo_dub": solo_dub,
            # "scrublet_pred": scrublet_pred,
            "dub_find_pred": dub_find_pred,
        }
        doublet_rates.append(row)
    dub_df = pd.DataFrame(doublet_rates)
    dub_df = dub_df.set_index("label")
    return dub_df


def lane_doublet_method_agreement(ad, root_path=None, est_per_1k=0.008):
    pass


# TODO: split out plot function
def remove_doublets(
    ad, root_path=None, solo_dub=False, solo_pred=False, dub_find=False, remove_scrublet=False,
):
    if not set(ALL_COLS).intersection(set(ad.obs.columns)):
        ad = _add_dub_methods(ad, root_path)
    # TODO: options to remove using others besides solo
    if dub_find:
        ad = ad[~(ad.obs["dub_find_pred"] == "True")]
    if solo_dub:
        ad = ad[~(ad.obs["solo_dub"] == "True")]
    if solo_pred:
        ad = ad[~(ad.obs["solo_pred"] == "True")]
    # if remove_scrublet:
    #     ad = ad[~ad.obs["solo_pred"]]
    return ad


def plot_doublets(ad, root_path=None):
    if not set(ALL_COLS).intersection(set(ad.obs.columns)):
        ad = _add_dub_methods(ad, root_path)
    cols = list(set(["n_genes", "total_counts"] + ALL_COLS).intersection(ad.obs.columns))
    sc.pl.umap(ad, color=cols)
