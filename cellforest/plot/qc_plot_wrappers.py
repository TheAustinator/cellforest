import matplotlib.pyplot as plt

from cellforest import CellBranch

DEFAULT_PLOT_RESOLUTION = (500, 500)  # width, height in pixels
PLOT_FILE_EXT = ".png"


def _get_plot_file_path(plot_func):
    plot_func_name = plot_func.__name__  # e.g., "plot_genes_per_cell_hist"
    plot_name = plot_func_name.replace("plot_", "", 1)  # e.g., "genes_per_cell_hist"

    plot_file_path = branch[run_name].plots_path / (
        plot_name + PLOT_FILE_EXT
    )  # inferred file path, e.g., "tests/data/example_usage/root/_plots/genes_per_cell_hist.png"

    # check if config requests overriding plot file name
    file_name_overridden = False
    plot_map = branch[run_name].plot_map
    for plot_map_key in [plot_name, plot_func_name]:
        if plot_map_key in plot_map and plot_map[plot_map_key].endswith(PLOT_FILE_EXT):
            plot_file_path = plot_map[plot_map_key]
            file_name_overridden = True

    if file_name_overridden:
        # TODO: properly log this information
        print(f"Using config-provided plot file name: {plot_file_path}")
    else:
        # TODO: properly log this information
        print(f"No viable plot file name was provided. Using default file name: {plot_file_path}")


def qc_plot_py(plot_func):
    def wrapper(branch: "CellBranch", **kwargs):
        fig = plt.gcf()
        dpi = fig.get_dpi()
        fig.set_size_inches(
            DEFAULT_PLOT_RESOLUTION[0] / float(dpi), DEFAULT_PLOT_RESOLUTION[1] / float(dpi)
        )  # scale to pixel resolution, irrespective of screen resolution

        plt.figure(plt.subplots(1, 1, figsize=(3, 3)))
        plot_func(branch, **kwargs)

        save_path = _get_plot_file_path(plot_func)
        plt.savefig(save_path, dpi=300)

    return wrapper


def qc_plot_r(plot_func):
    def wrapper(branch: "CellBranch", **kwargs):
        plot_func(branch, **kwargs)

    return wrapper
