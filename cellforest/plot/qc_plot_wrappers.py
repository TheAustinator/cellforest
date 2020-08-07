from functools import wraps
import matplotlib.pyplot as plt

from cellforest import CellBranch

DEFAULT_PLOT_RESOLUTION = (500, 500)  # width, height in pixels
PLOT_FILE_EXT = ".png"


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
        print(f"No viable plot file name was provided for {plot_name}, using default file path: {plot_file_path}")

    return plot_file_path


def qc_plot_py(plot_func):
    @wraps(plot_func)
    def wrapper(branch: "CellBranch", **kwargs):
        fig, ax = plt.subplots(1, 1)
        dpi = fig.get_dpi()
        fig.set_size_inches(
            DEFAULT_PLOT_RESOLUTION[0] / float(dpi), DEFAULT_PLOT_RESOLUTION[1] / float(dpi)
        )  # scale to pixel resolution, irrespective of screen resolution

        plot_func(branch, ax=ax, **kwargs)
        save_path = _get_plot_file_path(branch, plot_func)
        fig.savefig(save_path)

    return wrapper


def qc_plot_r(plot_func):
    @wraps(plot_func)
    def wrapper(branch: "CellBranch", **kwargs):
        plot_func(branch, **kwargs)

    return wrapper
