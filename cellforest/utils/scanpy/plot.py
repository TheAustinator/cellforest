from collections import deque
import math
from copy import deepcopy, copy
from multiprocessing import cpu_count
import logging
from numbers import Number as Num
from math import ceil
from pathlib import Path
from typing import Union, List, Optional, Tuple, Dict
import warnings

from IPython.core.display import display
from anndata import AnnData
from joblib import Parallel, delayed
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import spmatrix
import seaborn as sns
from matplotlib.axes import Axes
from seaborn.distributions import _DistributionPlotter

from cellforest.utils.scanpy.cell import marker_dict
from cellforest.utils.scanpy.gene import get_features, add_obsm

TEMPLATE_ATTRS_DEFAULT = (
    "collections",
    # "patches",
    # "lines",
    # "texts",
    # "artists",
    # "spines",
    # "xaxis",
    # "yaxis",
    # "tables",
    # "images",
    # "child_axes",
    # "legend_",
    # "patch",
)
# "title", "_left_title", "_right_title"


def _prep_axes(ax: Axes, rm_ticks: bool, rm_ax_labels: bool, tick_params: Optional[dict] = None):
    if rm_ticks:
        ax.set_xticks([])
        ax.set_yticks([])
    if rm_ax_labels:
        ax.set_xlabel("")
        ax.set_ylabel("")
    tick_params = tick_params if tick_params else dict()
    plt.tick_params(**tick_params)


def _get_kde_arr(feat: Union[np.ndarray, spmatrix], kde: Union[bool, np.ndarray, float]) -> np.ndarray:
    if isinstance(kde, bool):
        return np.array([kde] * feat.shape[1]).flatten()
    if isinstance(kde, float):
        fracs = np.array((feat > 0).sum(axis=0)).flatten()
        return fracs < kde * len(feat)
    if isinstance(kde, np.ndarray):
        return kde.flatten()


def _ax_copy_template(ax, ax_template, attrs: Tuple[str] = TEMPLATE_ATTRS_DEFAULT):
    for attr in attrs:
        try:
            setattr(ax, attr, deepcopy(getattr(ax_template, attr)))
        except Exception as e:
            if isinstance(getattr(ax_template, attr), (tuple, list, np.ndarray)):
                try:
                    setattr(ax, attr, list(map(copy, getattr(ax_template, attr))))
                except Exception as e:
                    setattr(ax, attr, copy(getattr(ax_template, attr)))
            else:
                setattr(ax, attr, copy(getattr(ax_template, attr)))


def _force_iter(x):
    if not isinstance(x, (list, set, tuple, np.ndarray, pd.Series)):
        x = [
            x,
        ]
    return x


def _replace_none(x, type_):
    return type_() if x is None else x


def _ax_copy(ax1):
    import pickle
    import io

    buf = io.BytesIO()
    pickle.dump(ax1, buf)
    buf.seek(0)
    ax2 = pickle.load(buf)
    return ax2


def _get_fig(
    n,
    titles: Optional[List[str]] = None,
    n_cols: int = 4,
    figsize: Optional[Tuple[float, float]] = None,
    flatten: bool = True,
    rm_ticks: bool = True,
    rm_ax_labels: bool = True,
    tick_params: Optional[dict] = None,
):
    titles = [titles,] if isinstance(titles, str) else titles
    n_cols = min(n, n_cols)
    n_rows = math.ceil(n / n_cols)
    ax_size = figsize if figsize else plt.rcParams["figure.figsize"]
    figsize = (ax_size * np.array([n_cols, n_rows])).tolist()
    fig, ax_arr = plt.subplots(n_rows, n_cols, figsize=figsize)
    if isinstance(ax_arr, np.ndarray) and flatten:
        ax_arr = ax_arr.flatten()
    elif isinstance(ax_arr, Axes):
        ax_arr = np.array([ax_arr])
    tick_prep_gen = map(lambda ax: _prep_axes(ax, rm_ticks, rm_ax_labels, tick_params), ax_arr)
    deque(tick_prep_gen, maxlen=0)
    if titles is not None:
        title_gen = map(lambda i, ax: ax.set_title(titles[i]), *zip(*enumerate(ax_arr[: len(titles)])))
        deque(title_gen, maxlen=0)
    return fig, ax_arr


def plot_gene_specificity(
    ad_list,
    genes,
    titles=None,
    n_cols=6,
    save_dir=None,
    save_prefix="",
    show=True,
    parallel=False,
    ignore_error=True,
    **kwargs,
):
    def _remove_ax_extras(ax):
        """remove asis labels and ticks"""
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_xlabel(None)
        ax.set_ylabel(None)

    def _single_gene_plot(g):
        n_rows = math.ceil(n_plots / n_cols)
        fig, ax_list = plt.subplots(n_rows, n_cols, figsize=(3.3 * n_cols, 3 * n_rows))
        plt.tick_params(
            axis="x",  # changes apply to the x-axis
            which="both",  # both major and minor ticks are affected
            bottom=False,  # ticks along the bottom edge are off
            top=False,  # ticks along the top edge are off
            labelbottom=False,  # labels along the bottom edge are off
        )
        if isinstance(ax_list, np.ndarray):
            ax_list = ax_list.flatten()
        for i, ad in enumerate(ad_list):
            get_kwarg = lambda v: v[i] if isinstance(v, list) else v
            _kwargs = {k: get_kwarg(v) for k, v in kwargs.items()}
            title = titles[i] if titles is not None else titles
            i *= 2
            (ax1, ax2) = ax_list[i : i + 2]
            _remove_ax_extras(ax1)
            _remove_ax_extras(ax2)
            umap_kde_scatter(ad, title=title, var=g, ax_scat=ax1, ax_dens=ax2, show=False, **_kwargs)
        if save_dir:
            path = Path(save_dir) / f"{save_prefix}{g}.png"
            path.parent.mkdir(parents=True, exist_ok=True)
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
            try:
                _single_gene_plot(g)
            except Exception as e:
                logging.warning(
                    f"Error encountered for gene: {g}. Error will be raised if `ignore_error`. You have: `ignore_error`={ignore_error}"
                )
                if not ignore_error:
                    raise e


def umap(
    adata,
    color=None,
    color_obsm=None,
    n_obsm=None,
    group: Optional[str] = None,
    vmin_frac=None,
    vmax_frac=0.9,
    vmin_frac_nonzero=False,
    vmax_frac_nonzero=True,
    raw=True,
    layer="norm",
    kde: Union[bool, float] = 0.1,
    ax: Optional[Axes] = None,
    show: bool = False,
    return_fig: bool = False,
    ignore_missing: bool = False,
    **kwargs,
):
    """
    Wrapper for scanpy umap plot with additional functionality
    Args:
        adata:
        color:
        color_obsm:
        n_obsm:
        group: obs column by which to stratify kde lines in kde plots
        vmin_frac: fractional percentile at which to set colorbar min
        vmax_frac: fractional percentile at which to sedef umapt colorbar max
        vmin_frac_nonzero:
        vmax_frac_nonzero:
        raw:
        kde: If boolean, whether to make all plots kde rather than scatter. If
            float, the minimum fraction of expressing cells that triggers a
            switch from scatter to kde for better visibility
        kwargs:
    """
    feat = None
    color = list() if color is None else color
    color = [color,] if isinstance(color, str) else color
    color = list(color)
    if ignore_missing:
        missing = set(color).difference(adata.var_names, adata.obs_names)
        logging.warning(f"`ignore_missing=True`, ignoring {len(missing)}/{len(color)}")
        color = list(set(color).intersection(adata.var_names, adata.obs_names))

    def _get_v(v_frac, _feat, nonzero):
        """get vmin or vmax based on percentile"""
        if isinstance(v_frac, Num):
            v_frac = np.repeat(v_frac, len(color))
        v_frac = np.array(v_frac)
        _feat = _feat
        if nonzero:
            frac_zero = (_feat == 0).sum(axis=0) / _feat.shape[0]
            v_frac = frac_zero + (1 - frac_zero) * v_frac
        v = np.diagonal(np.percentile(_feat, v_frac * 100, axis=0)).tolist()
        return v

    if raw or layer:
        adata = adata.copy()
    if raw:
        if adata.raw is None:
            logging.warning(
                "Attempted to use `raw` attribute for umap, but `raw` attribute missing. Set `raw=False` to silence."
            )
        adata.X = adata.raw.X
    if layer:
        if layer in adata.layers:
            adata.X = adata.layers[layer]
        else:
            default = "raw" if raw else "X"
            logging.warning(f"layer: {layer} not present in layers. Defaulting to {default}")
    if (color or color_obsm) and (vmin_frac or vmax_frac):
        color_obsm = _force_iter(color_obsm) if color_obsm else []
        for obsm_key in color_obsm:
            color += adata.obs.columns[adata.obs.columns.str.startswith(obsm_key)].tolist()
        adata = add_obsm(adata, color_obsm + ["X_umap"], n_obsm, copy=True)
        feat = get_features(adata, color, color_obsm, n_obsm, ignore_missing=ignore_missing)
        if vmin_frac is not None:
            kwargs["vmin"] = _get_v(vmin_frac, feat, vmin_frac_nonzero)
        if vmax_frac is not None:
            kwargs["vmax"] = _get_v(vmax_frac, feat, vmax_frac_nonzero)
    if kde is False:
        return sc.pl.umap(adata, color=color, **kwargs)
    feat = feat if feat is not None else get_features(adata, color, color_obsm, n_obsm, ignore_missing=ignore_missing)
    kde_arr = _get_kde_arr(feat, kde)
    if not kde_arr.any():
        return None, sc.pl.umap(adata, color=color, ax=ax, **kwargs)
    if not len(color) == len(kde_arr):
        print(kde_arr.flatten().shape)
        raise ValueError(f"Length mismatch between `key` ({len(color)}) and `kde_arr` ({len(kde_arr)})")
    feat = feat.toarray() if isinstance(feat, spmatrix) else feat
    df = pd.DataFrame(feat, columns=color)
    obs_names = ["X_umap_0", "X_umap_1"]
    df_umap = get_features(adata, obs_names, f_obsm="X_umap", ignore_missing=ignore_missing)
    df_grp = get_features(adata, group, ignore_missing=ignore_missing)
    df = pd.concat([df, df_umap, df_grp], axis=1)
    if not ax:
        fig, ax_arr = _get_fig(len(color), color,)
    else:
        fig, ax_arr = None, []
    # ax_template = None
    for i, x in enumerate(color):
        _ax = ax if ax else ax_arr[i]
        if not kde_arr[i]:
            _kwargs = copy(kwargs)
            for v_name in ("vmin", "vmax"):
                v_val = _kwargs.pop(v_name, None)
                if isinstance(v_val, (list, tuple, np.ndarray)):
                    v_val = v_val[i]
                _kwargs[v_name] = v_val
            sc.pl.umap(adata, color=x, ax=_ax, **_kwargs, show=False)
        else:
            # if not ax_template:
            #     ax_template = copy(_single_kde(df, group=group, ax=ax, **kwargs))
            # else:
            #     _ax_copy_template(ax, ax_template)
            # ax.update_from(ax_template)
            _single_kde(df, group=group, ax=_ax, **kwargs)
            _single_sparse_scatter(df, color=color[i], ax=_ax, **kwargs)
            _ax.set_title(kwargs.get("title", color[i]))
    if show:
        plt.show()
    if return_fig:
        return fig, ax_arr


def _single_kde(
    df,
    group: Optional[str] = None,
    kde_levels=7,
    kde_alpha=0.5,
    figsize=(12, 5),
    ax: Optional[Axes] = None,
    kde_kwargs=None,
    **kwargs,
):

    kde_kwargs = _replace_none(kde_kwargs, dict)
    if not ax:
        figsize = figsize if figsize else plt.rcParams["figure.figsize"]
        fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax = sns.kdeplot(
        data=df,
        x="X_umap_0",
        y="X_umap_1",
        hue=group,
        alpha=kde_alpha,
        levels=kde_levels,
        legend=True,
        ax=ax,
        **kde_kwargs,
    )
    return ax


def _single_sparse_scatter(
    df, color=None, cmap="flare", scat_alpha=0.5, dot_scale=1, ax: Optional[Axes] = None, scat_kwargs=None, **kwargs
):
    scat_kwargs = _replace_none(scat_kwargs, dict)
    df_g = df[df[color] > 0]
    norm = plt.Normalize(0, df_g[color].max())
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    s = dot_scale * 10 / np.log(len(df_g) / 1000 + 2)
    if isinstance(cmap, str):
        cmap = sns.color_palette(cmap, as_cmap=True)
    sns.scatterplot(
        x="X_umap_0",
        y="X_umap_1",
        hue=color,
        s=s,
        data=df_g,
        cmap=cmap,
        alpha=scat_alpha,
        ax=ax,
        legend=False,
        **scat_kwargs,
    )
    return ax


def kdeplot(
    x=None,  # Allow positional x, because behavior will not change with reorg
    *,
    y=None,
    shade=None,  # Note "soft" deprecation, explained below
    vertical=False,  # Deprecated
    kernel=None,  # Deprecated
    bw=None,  # Deprecated
    gridsize=200,  # TODO maybe depend on uni/bivariate?
    cut=3,
    clip=None,
    legend=True,
    cumulative=False,
    shade_lowest=None,  # Deprecated, controlled with levels now
    cbar=False,
    cbar_ax=None,
    cbar_kws=None,
    ax=None,
    # New params
    weights=None,  # TODO note that weights is grouped with semantics
    hue=None,
    palette=None,
    hue_order=None,
    hue_norm=None,
    multiple="layer",
    common_norm=True,
    common_grid=False,
    levels=10,
    thresh=0.05,
    bw_method="scott",
    bw_adjust=1,
    log_scale=None,
    color=None,
    fill=None,
    # Renamed params
    data=None,
    data2=None,
    **kwargs,
):
    """FROM SEABORN"""
    # Handle deprecation of `data2` as name for y variable
    if data2 is not None:

        y = data2

        # If `data2` is present, we need to check for the `data` kwarg being
        # used to pass a vector for `x`. We'll reassign the vectors and warn.
        # We need this check because just passing a vector to `data` is now
        # technically valid.

        x_passed_as_data = x is None and data is not None and np.ndim(data) == 1

        if x_passed_as_data:
            msg = "Use `x` and `y` rather than `data` `and `data2`"
            x = data
        else:
            msg = "The `data2` param is now named `y`; please update your code"

        warnings.warn(msg, FutureWarning)

    # Handle deprecation of `vertical`
    if vertical:
        msg = (
            "The `vertical` parameter is deprecated and will be removed in a "
            "future version. Assign the data to the `y` variable instead."
        )
        warnings.warn(msg, FutureWarning)
        x, y = y, x

    # Handle deprecation of `bw`
    if bw is not None:
        msg = (
            "The `bw` parameter is deprecated in favor of `bw_method` and "
            f"`bw_adjust`. Using {bw} for `bw_method`, but please "
            "see the docs for the new parameters and update your code."
        )
        warnings.warn(msg, FutureWarning)
        bw_method = bw

    # Handle deprecation of `kernel`
    if kernel is not None:
        msg = "Support for alternate kernels has been removed. " "Using Gaussian kernel."
        warnings.warn(msg, UserWarning)

    # Handle deprecation of shade_lowest
    if shade_lowest is not None:
        if shade_lowest:
            thresh = 0
        msg = (
            "`shade_lowest` is now deprecated in favor of `thresh`. "
            f"Setting `thresh={thresh}`, but please update your code."
        )
        warnings.warn(msg, UserWarning)

    # Handle `n_levels`
    # This was never in the formal API but it was processed, and appeared in an
    # example. We can treat as an alias for `levels` now and deprecate later.
    levels = kwargs.pop("n_levels", levels)

    # Handle "soft" deprecation of shade `shade` is not really the right
    # terminology here, but unlike some of the other deprecated parameters it
    # is probably very commonly used and much hard to remove. This is therefore
    # going to be a longer process where, first, `fill` will be introduced and
    # be used throughout the documentation. In 0.12, when kwarg-only
    # enforcement hits, we can remove the shade/shade_lowest out of the
    # function signature all together and pull them out of the kwargs. Then we
    # can actually fire a FutureWarning, and eventually remove.
    if shade is not None:
        fill = shade

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

    p = _DistributionPlotter(data=data, variables=_DistributionPlotter.get_semantics(locals()),)

    p.map_hue(palette=palette, order=hue_order, norm=hue_norm)

    if ax is None:
        ax = plt.gca()

    p._attach(ax, allowed_types=["numeric", "datetime"], log_scale=log_scale)

    method = ax.fill_between if fill else ax.plot

    color = _default_color(method, hue, color, kwargs)

    if not p.has_xy_data:
        return ax

    # Pack the kwargs for statistics.KDE
    estimate_kws = dict(
        bw_method=bw_method, bw_adjust=bw_adjust, gridsize=gridsize, cut=cut, clip=clip, cumulative=cumulative,
    )
    p.plot_bivariate_density(
        common_norm=common_norm,
        fill=fill,
        levels=levels,
        thresh=thresh,
        legend=legend,
        color=color,
        cbar=cbar,
        cbar_ax=cbar_ax,
        cbar_kws=cbar_kws,
        estimate_kws=estimate_kws,
        **kwargs,
    )


def _default_color(method, hue, color, kws):
    """FROM SEABORN. If needed, get a default key by using the matplotlib property cycle."""
    from matplotlib.colors import to_rgb

    if hue is not None:
        # This warning is probably user-friendly, but it's currently triggered
        # in a FacetGrid context and I don't want to mess with that logic right now
        #  if key is not None:
        #      msg = "`key` is ignored when `hue` is assigned."
        #      warnings.warn(msg)
        return None

    if color is not None:
        return color

    elif method.__name__ == "plot":

        (scout,) = method([], [], **kws)
        color = scout.get_color()
        scout.remove()

    elif method.__name__ == "scatter":

        # Matplotlib will raise if the size of x/y don't match s/c,
        # and the latter might be in the kws dict
        scout_size = max(
            np.atleast_1d(kws.get(key, [])).shape[0] for key in ["s", "c", "fc", "facecolor", "facecolors"]
        )
        scout_x = scout_y = np.full(scout_size, np.nan)

        scout = method(scout_x, scout_y, **kws)
        facecolors = scout.get_facecolors()

        if not len(facecolors):
            # Handle bug in matplotlib <= 3.2 (I think)
            # This will limit the ability to use non key= kwargs to specify
            # a key in versions of matplotlib with the bug, but trying to
            # work out what the user wanted by re-implementing the broken logic
            # of inspecting the kwargs is probably too brittle.
            single_color = False
        else:
            single_color = np.unique(facecolors, axis=0).shape[0] == 1

        # Allow the user to specify an array of colors through various kwargs
        if "c" not in kws and single_color:
            color = to_rgb(facecolors[0])

        scout.remove()

    elif method.__name__ == "bar":

        # bar() needs masked, not empty data, to generate a patch
        (scout,) = method([np.nan], [np.nan], **kws)
        color = to_rgb(scout.get_facecolor())
        scout.remove()

    elif method.__name__ == "fill_between":

        # There is a bug on matplotlib < 3.3 where fill_between with
        # datetime units and empty data will set incorrect autoscale limits
        # To workaround it, we'll always return the first key in the cycle.
        # https://github.com/matplotlib/matplotlib/issues/17586
        ax = method.__self__
        datetime_axis = any(
            [
                isinstance(ax.xaxis.converter, mpl.dates.DateConverter),
                isinstance(ax.yaxis.converter, mpl.dates.DateConverter),
            ]
        )

        kws = _normalize_kwargs(kws, mpl.collections.PolyCollection)

        scout = method([], [], **kws)
        facecolor = scout.get_facecolor()
        color = to_rgb(facecolor[0])
        scout.remove()

    return color


def _normalize_kwargs(kws, artist):
    """FROM SEABORN. Wrapper for mpl.cbook.normalize_kwargs that supports <= 3.2.1."""

    from matplotlib.cbook import normalize_kwargs

    _alias_map = {
        "key": ["c"],
        "linewidth": ["lw"],
        "linestyle": ["ls"],
        "facecolor": ["fc"],
        "edgecolor": ["ec"],
        "markerfacecolor": ["mfc"],
        "markeredgecolor": ["mec"],
        "markeredgewidth": ["mew"],
        "markersize": ["ms"],
    }
    try:
        kws = normalize_kwargs(kws, artist)
    except AttributeError:
        kws = normalize_kwargs(kws, _alias_map)
    return kws


# def _umap_kde(
#     ad: AnnData, key: Optional[Union[str, List[str]]], n_cols: int = 4, color_obsm=None, n_obsm=None,
# ):
#     if isinstance(key, str):
#         key = [
#             key,
#         ]
#     n_rows = ceil(len(key) / n_cols)
#     figsize = (np.array(plt.rcParams["figure.figsize"]) * [n_rows, n_cols]).tolist()
#     fig, ax_arr = plt.subplots(n_rows, n_cols, figsize=figsize)
#     if isinstance(ax_arr, np.ndarray):
#         ax_arr = ax_arr.flatten()
#     feat = get_features(ad, key, color_obsm, n_obsm)
#     df = pd.DataFrame(feat, columns=key)
#
#     sns.kdeplot(
#         data=df,
#         x="UMAP1",
#         y="UMAP2",
#         hue=kde_hue,
#         alpha=kde_alpha,
#         levels=kde_levels,
#         legend=True,
#         ax=ax,
#         **kde_kwargs,
#     )


def umap_kde_scatter(
    ad,
    var=(),
    obs=(),
    groupby=None,
    title=None,
    vmin_frac=0,
    vmax_frac=0.95,
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
    """DEPRECATE: SUNSET THIS"""

    def replace_none(x, type_):
        return type_() if x is None else x

    if kde_hue not in obs:
        obs = list(obs) + [kde_hue]
    kde_kwargs = replace_none(kde_kwargs, dict)
    scat_kwargs = replace_none(scat_kwargs, dict)
    df = _get_ad_df(ad, var, obs)

    for g in var:
        if not (ax_scat and ax_dens):
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
        else:
            (ax1, ax2) = (ax_scat, ax_dens)
        sns.kdeplot(
            data=df,
            x="UMAP1",
            y="UMAP2",
            hue=kde_hue,
            alpha=kde_alpha,
            levels=kde_levels,
            legend=True,
            ax=ax2,
            **kde_kwargs,
        )
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
        umap(ad, color=g, vmin_frac=vmin_frac, vmax_frac=vmax_frac, ax=ax1, show=show, **scat_kwargs)
        title = f"{title} {g}" if title is not None else g
        ax1.set_title(title)
        ax2.set_title(title)
        if show:
            plt.show()


def _get_ad_df(ad, var, obs):
    var = list(_force_iter(var))
    obs = list(_force_iter(obs))
    umap_arr = ad.obsm["X_umap"]
    var_arr = ad[:, var].X.toarray()
    obs_arr = ad.obs[obs].to_numpy()
    arr = [umap_arr, var_arr, obs_arr]
    arr = np.concatenate([_arr for _arr in arr if np.multiply(*_arr.shape)], axis=1)
    cols = ("UMAP1", "UMAP2", *var, *obs)
    df = pd.DataFrame(arr, columns=cols)
    df[["UMAP1", "UMAP2"]] = df[["UMAP1", "UMAP2"]].astype(float)
    return df


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


def marker_matrix(
    ad: AnnData,
    gt_thresh: Optional[Dict[str, float]] = None,
    lt_thresh: Optional[Dict[str, float]] = None,
    sort_by: str = "volcano",
    sort_ascending: bool = False,
    max_n: Optional[int] = None,
    prefix_filter: Optional[Union[List[str], str]] = None,
    groupby: Optional[str] = None,
    use_raw: bool = False,
    layer: Optional[str] = "norm",
    dendrogram: bool = True,
    standard_scale: str = "var",
    **kwargs,
):
    genes_dict, _ = marker_dict(
        ad,
        gt_thresh=gt_thresh,
        lt_thresh=lt_thresh,
        sort_by=sort_by,
        sort_ascending=sort_ascending,
        max_n=max_n,
        prefix_filter=prefix_filter,
    )
    groupby = groupby if groupby else ad.uns["rank_genes_groups"]["params"]["groupby"]
    sc.pl.matrixplot(
        ad,
        var_names=genes_dict,
        groupby=groupby,
        use_raw=use_raw,
        layer=layer,
        dendrogram=dendrogram,
        standard_scale=standard_scale,
        **kwargs,
    )
