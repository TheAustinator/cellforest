from functools import wraps
import json
from pathlib import Path
import logging

from IPython.display import Image
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cellforest import CellBranch

_DEFAULT_PLOT_RESOLUTION_PX = (500, 500)  # width, height in pixels
_DEFAULT_BIG_PLOT_RESOLUTION_PX = (1000, 1000)  # width, height in pixels
_PLOT_FILE_EXT = ".png"

_R_PLOT_SCRIPTS_PATH = Path(__file__).parent / "r"
_R_FUNCTIONS_FILEPATH = Path(__file__).parent.parent / "processes/scripts/functions.R"

_NONE_VARIATIONS = [None, "none", "None", "NULL", "NA"]

_LOG = logging.getLogger("qc_plot_wrappers")


# TODO: move these bad boiz to dataforest
def qc_plot_py(plot_func):
    @wraps(plot_func)
    def wrapper(branch: "CellBranch", **kwargs):
        plot_path = kwargs.pop("plot_path", None)
        stratify = kwargs.pop("stratify", None)
        fig, ax = _prep_plot(kwargs)

        if plot_path is not None:
            matplotlib.use("Agg")  # don't plot on screen
        if stratify not in _NONE_VARIATIONS:
            try:
                labels = branch.meta[stratify]
                if isinstance(labels, pd.DataFrame):
                    cols = list(labels)
                    labels = labels[cols].apply(lambda row: "_".join(row.values.astype(str)), axis=1)
                    stratify = "_".join(stratify)
                    branch = branch.copy()
                    branch.meta[stratify] = labels
                kwargs["labels"] = labels
                kwargs["legend_title"] = stratify
            except KeyError:
                _LOG.warning(f"{plot_func.__name__} with key '{stratify}' is skipped because key is not in metadata")
                return
        try:
            # "labels" key only works for `Counts` methods
            plot_func(branch, ax=ax, **kwargs)
        except AttributeError as e:
            if stratify in _NONE_VARIATIONS:
                raise e
            ax.clear()
            _stratified_plot(plot_func, stratify, branch, ax, **kwargs)
        if plot_path is not None:
            _LOG.info(f"saving py figure to {plot_path}")
            fig.savefig(plot_path)
        return fig, ax

    return wrapper


def qc_plot_r(plot_func):
    @wraps(plot_func)
    def wrapper(branch: "CellBranch", **kwargs):
        # TODO: move temp spec to a hook
        r_script = _R_PLOT_SCRIPTS_PATH / (plot_func.__name__ + ".R")
        plot_size = kwargs.pop("plot_size", _DEFAULT_PLOT_RESOLUTION_PX)
        stratify = kwargs.pop("stratify", None)
        plot_path = kwargs.pop("plot_path", None)

        if stratify not in _NONE_VARIATIONS:
            if stratify in branch.meta:  # column exists in metadata
                kwargs["group.by"] = stratify
            else:
                _LOG.warning(f"{plot_func.__name__} with key '{stratify}' is skipped because key is not in metadata")
                return

        args = [  # corresponding arguments in r/plot_entry_point.R
            _R_PLOT_SCRIPTS_PATH,  # r_plot_scripts_path
            branch.paths["root"],  # root_dir
            branch.spec.shell_str,  # spec_str
            branch.current_process,  # current_process
            plot_path,  # plot_file_path
            plot_size[0],  # plot_width_px
            plot_size[1],  # plot_height_px
            _R_FUNCTIONS_FILEPATH,  # r_functions_filepath
            json.dumps("kwargs = " + str(kwargs if kwargs else {})),  # TODO-QC: is there a better way to handle this?
        ]
        _LOG.info(f"saved R figure to {plot_path}")
        plot_func(branch, r_script, args)  # kwargs already included in args
        return Image("/tmp/plot.png")

    return wrapper


# noinspection PyPep8Naming
class requires:
    def __init__(self, req_process):
        self._req_process = req_process

    def __call__(self, func):
        @wraps(func)
        def wrapper(branch: "CellBranch", *args, **kwargs):
            precursors = branch.spec.get_precursors_lookup(incl_current=True)[branch.current_process]
            if self._req_process not in precursors:
                proc = self._req_process
                raise ValueError(
                    f"This plot method requires a branch at `{proc}` or later. Current process run: {precursors}. If "
                    f"`{proc}` has already been run, please use `branch.goto_process`. Otherwise, please run `{proc}`."
                )
            return func(branch, *args, **kwargs)

        return wrapper


def _prep_plot(kwargs):
    xlim = kwargs.pop("xlim", None)
    ylim = kwargs.pop("ylim", None)
    xscale = kwargs.pop("xscale", None)
    yscale = kwargs.pop("yscale", None)
    plot_size = kwargs.pop("plot_size", _DEFAULT_PLOT_RESOLUTION_PX)
    figsize = kwargs.pop("figsize", None)
    if "ax" in kwargs:
        ax = kwargs.pop("ax")
        fig = kwargs.pop("fig", plt.gcf())
    else:
        fig, ax = plt.subplots(1, 1)
    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)
    if xscale:
        ax.set_xscale(xscale)
    if yscale:
        ax.set_yscale(yscale)
    dpi = fig.get_dpi()
    # scale to pixel resolution, irrespective of screen resolution
    fig.set_size_inches(plot_size[0] / float(dpi), plot_size[1] / float(dpi))
    if figsize:
        fig.set_size_inches(*figsize)
    return fig, ax


def _stratified_plot(plot_func, stratify, branch, ax, **kwargs):
    """
    labels = np.unique(np.array(kwargs.pop("labels")))
    legend_title = kwargs.pop("legend_title", None)
    ax.legend(title=legend_title)
    for label in labels:
        branch_sub = branch.copy()
        meta_sub = branch.meta[branch.meta[stratify] == label]
        branch_sub.set_meta(meta_sub)
        plot_func(branch_sub, ax=ax, label=label, **kwargs)
    """
    labels = np.unique(np.array(kwargs.pop("labels")))
    legend_title = kwargs.pop("legend_title", None)
    for label in labels:
        branch_sub = branch.copy()
        meta_sub = branch.meta[branch.meta[stratify] == label]
        branch_sub.set_meta(meta_sub)
        plot_func(branch_sub, label=label, **kwargs)
        ax.legend(title=legend_title)
