from os import remove
import pickle
from functools import wraps
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt

from cellforest import CellBranch

matplotlib.use("Agg")

DEFAULT_PLOT_RESOLUTION_PX = (500, 500)  # width, height in pixels
PLOT_FILE_EXT = ".png"

R_PLOT_SCRIPTS_PATH = Path(__file__).parent / "r"
R_FUNCTIONS_FILEPATH = Path(__file__).parent.parent / "processes/scripts/functions.R"


def _get_plot_file_path(branch: "CellBranch", plot_func):
    plot_func_name = plot_func.__name__  # e.g., "plot_genes_per_cell_hist"
    plot_name = plot_func_name.replace("plot_", "", 1)  # e.g., "genes_per_cell_hist"

    run_name = branch.current_process
    plot_file_path = branch[run_name].plots_path / (
        plot_name + PLOT_FILE_EXT
    )  # default file path, e.g., "tests/data/example_usage/root/_plots/genes_per_cell_hist.png"

    # check if config requests overriding plot file name
    file_name_overridden = False
    plot_map = branch[run_name].plot_map
    for plot_map_key in [plot_name, plot_func_name]:
        if plot_map_key in plot_map and str(plot_map[plot_map_key]).endswith(PLOT_FILE_EXT):
            plot_file_path = plot_map[plot_map_key]
            file_name_overridden = True

    if file_name_overridden:
        # TODO: properly log this information
        print(f"Using provided plot file path in config: {plot_file_path}")
    else:
        # TODO: properly log this information
        print(f"No valid plot file path provided for {plot_name}, using default file path: {plot_file_path}")

    return plot_file_path


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
        fig, ax = plt.subplots(1, 1)
        dpi = fig.get_dpi()
        fig.set_size_inches(
            DEFAULT_PLOT_RESOLUTION_PX[0] / float(dpi), DEFAULT_PLOT_RESOLUTION_PX[1] / float(dpi)
        )  # scale to pixel resolution, irrespective of screen resolution

        plot_func(branch, ax=ax, **kwargs)
        save_path = _get_plot_file_path(branch, plot_func)
        fig.savefig(save_path)

    return wrapper


def qc_plot_r(plot_func):
    @wraps(plot_func)
    def wrapper(branch: "CellBranch", **kwargs):
        # TODO: move temp spec to a hook
        temp_spec_path = _create_temp_spec(branch)
        r_script = R_PLOT_SCRIPTS_PATH / (plot_func.__name__ + ".R")
        save_path = _get_plot_file_path(branch, plot_func)

        args = [  # corresponding arguments in r/plot_entry_point.R
            branch.paths["root"],  # root_dir
            temp_spec_path,  # path_to_temp_spec
            branch.current_process,  # current_process
            save_path,  # plot_file_path
            R_PLOT_SCRIPTS_PATH,  # r_plot_scripts_path
            DEFAULT_PLOT_RESOLUTION_PX[0],  # plot_width_px
            DEFAULT_PLOT_RESOLUTION_PX[1],  # plot_height_px
            R_FUNCTIONS_FILEPATH,  # r_functions_filepath
        ]

        plot_func(branch, r_script, args, **kwargs)
        _remove_temp_spec(temp_spec_path)

    return wrapper
