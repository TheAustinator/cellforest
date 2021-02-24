from multiprocessing import cpu_count
from pathlib import Path

from IPython.core.display import display
from anndata import AnnData
from joblib import Parallel, delayed
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns


def plot_gene_specificity(ad_list, genes, save_dir=None, save_prefix="", show=True, parallel=False, **kwargs):
    def _single_gene_plot(g):
        path = Path(save_dir) / f"{save_prefix}{g}.png"
        path.parent.mkdir(parents=True, exist_ok=True)
        fig, ax_list = plt.subplots(1, n_plots, figsize=(3.3 * n_plots, 3))
        for i, ad in enumerate(ad_list):
            get_kwarg = lambda v: v[i] if isinstance(v, list) else v
            _kwargs = {k: get_kwarg(v) for k, v in kwargs.items()}
            i *= 2
            (ax1, ax2) = ax_list[i : i + 2]
            umap_kde(ad, var=g, ax_scat=ax1, ax_dens=ax2, show=False, **_kwargs)
        if save_dir:
            fig.savefig(path)
        if show:
            display(fig)

    if isinstance(genes, str):
        genes = (genes,)
    n_plots = 2 * len(ad_list)
    if parallel:
        raise NotImplemented("AnnData currently breaks parallelism")
        Parallel(n_jobs=cpu_count() - 1)(delayed(_single_gene_plot)(g) for g in genes)
    else:
        for g in genes:
            _single_gene_plot(g)


def umap_kde(
    ad,
    var=(),
    obs=(),
    groupby=None,
    cmap="flare",
    kde_hue=None,
    kde_alpha=0.5,
    kde_levels=7,
    scat_alpha=0.5,
    dot_scale=1,
    figsize=(12, 5),
    show=True,
    ax_scat=None,
    ax_dens=None,
    kde_kwargs=None,
    scat_kwargs=None,
):
    def force_iter(x):
        if not isinstance(x, (list, set, tuple, np.ndarray, pd.Series)):
            x = (x,)
        return x

    def replace_none(x, type_):
        return type_() if x is None else x

    var = list(force_iter(var))
    obs = list(force_iter(obs))
    if kde_hue not in obs:
        obs.append(kde_hue)
    kde_kwargs = replace_none(kde_kwargs, dict)
    scat_kwargs = replace_none(scat_kwargs, dict)
    umap_arr = ad.obsm["X_umap"]
    var_arr = ad[:, var].X.toarray()
    obs_arr = ad.obs[obs].to_numpy()
    arr = [umap_arr, var_arr, obs_arr]
    arr = np.concatenate([_arr for _arr in arr if np.multiply(*_arr.shape)], axis=1)
    cols = ("UMAP1", "UMAP2", *var, *obs)
    df = pd.DataFrame(arr, columns=cols)
    df[["UMAP1", "UMAP2"]] = df[["UMAP1", "UMAP2"]].astype(float)
    for g in var:
        if not (ax_scat and ax_dens):
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        else:
            (ax1, ax2) = (ax_scat, ax_dens)
        sns.kdeplot(data=df, x="UMAP1", y="UMAP2", hue=kde_hue, alpha=kde_alpha, levels=kde_levels, ax=ax2)
        df_g = df[df[g] > 0]
        norm = plt.Normalize(0, df_g[g].max())
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm.set_array([])
        s = dot_scale * 10 / np.log(len(df_g) / 1000 + 2)
        if isinstance(cmap, str):
            cmap = sns.color_palette(cmap, as_cmap=True)
        sns.scatterplot(x="UMAP1", y="UMAP2", hue=g, s=s, data=df_g, cmap=cmap, alpha=scat_alpha, ax=ax2)
        ax2.get_legend().remove()
        if show:
            ax2.figure.colorbar(sm)
        ax2.set_title(g)
        sc.pl.umap(ad, color=g, ax=ax1, show=show)
        if show:
            plt.show()
