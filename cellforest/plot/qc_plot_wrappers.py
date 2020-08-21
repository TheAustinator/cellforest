from os import remove
import pickle
from functools import wraps
import json
from functools import update_wrapper
from pathlib import Path
import logging

import matplotlib
import matplotlib.pyplot as plt

from cellforest import CellBranch

DEFAULT_PLOT_RESOLUTION_PX = (500, 500)  # width, height in pixels
DEFAULT_BIG_PLOT_RESOLUTION_PX = (1000, 1000)  # width, height in pixels
PLOT_FILE_EXT = ".png"

R_PLOT_SCRIPTS_PATH = Path(__file__).parent / "r"
R_FUNCTIONS_FILEPATH = Path(__file__).parent.parent / "processes/scripts/functions.R"

DEFAULT_ASSAY = "rna"
NONE_VARIATIONS = [None, "none", "None", "NULL", "NA"]


def _create_temp_spec(branch: "CellBranch"):
    """
    Create a temporary text file with spec as a string to pass to R
    """
    run_name = branch.current_process
    plots_path = branch[run_name].plots_path
    temp_spec_path = plots_path / "temp_spec"

    with open(str(temp_spec_path), "wb") as f:
        pickle.dump(branch.spec, f, protocol=pickle.HIGHEST_PROTOCOL)

    return temp_spec_path


def _remove_temp_spec(temp_spec_path: str):
    remove(temp_spec_path)


def qc_plot_py(plot_func):
    @wraps(plot_func)
    def wrapper(branch: "CellBranch", **kwargs):
        matplotlib.use("Agg")  # don't plot on screen
        plot_size = kwargs.pop("plot_size", DEFAULT_PLOT_RESOLUTION_PX)
        stratify = kwargs.pop("stratify", None)
        plot_path = kwargs.pop("plot_path", None)

        fig, ax = plt.subplots(1, 1)
        dpi = fig.get_dpi()
        fig.set_size_inches(
            plot_size[0] / float(dpi), plot_size[1] / float(dpi)
        )  # scale to pixel resolution, irrespective of screen resolution

        if stratify not in NONE_VARIATIONS:
            try:
                kwargs["labels"] = branch.meta[stratify]
                kwargs["legend_title"] = stratify
            except KeyError:
                logging.warning(f"{plot_func.__name__} with key '{stratify}' is skipped because key is not in metadata")
                return
        else:
            kwargs["labels"] = [DEFAULT_ASSAY] * len(branch.meta)

        plot_func(branch, ax=ax, **kwargs)
        fig.savefig(plot_path)

    return wrapper


def qc_plot_r(plot_func):
    @wraps(plot_func)
    def wrapper(branch: "CellBranch", **kwargs):
        # TODO: move temp spec to a hook
        temp_spec_path = _create_temp_spec(branch)
        r_script = R_PLOT_SCRIPTS_PATH / (plot_func.__name__ + ".R")
        plot_size = kwargs.pop("plot_size", DEFAULT_PLOT_RESOLUTION_PX)
        stratify = kwargs.pop("stratify", None)
        plot_path = kwargs.pop("plot_path", None)

        if stratify not in NONE_VARIATIONS:
            if stratify in branch.meta:  # column exists in metadata
                kwargs["group.by"] = stratify
            else:
                logging.warning(f"{plot_func.__name__} with key '{stratify}' is skipped because key is not in metadata")
                return

        args = [  # corresponding arguments in r/plot_entry_point.R
            R_PLOT_SCRIPTS_PATH,  # r_plot_scripts_path
            branch.paths["root"],  # root_dir
            temp_spec_path,  # path_to_temp_spec
            branch.current_process,  # current_process
            plot_path,  # plot_file_path
            plot_size[0],  # plot_width_px
            plot_size[1],  # plot_height_px
            R_FUNCTIONS_FILEPATH,  # r_functions_filepath
            json.dumps("kwargs = " + str(kwargs if kwargs else {})),  # TODO-QC: is there a better way to handle this?
        ]

        plot_func(branch, r_script, args)  # kwargs already included in args
        _remove_temp_spec(temp_spec_path)

    return wrapper
