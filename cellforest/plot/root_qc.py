import matplotlib.pyplot as plt

from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_process_r_script

ALL_QC_PLOTS_ROOT = [  # reference of all supported QC plots
    "genes_per_cell_hist",
    "umis_per_cell_hist",
    "umis_vs_genes_hist",
    "perc_mito_per_cell_hist",
    "highest_exprs_genes_dens_plt",
]

DEFAULT_QC_PLOTS_ROOT = [  # default QC plots
    "genes_per_cell_hist",
    "umis_per_cell_hist",
    "umis_vs_genes_hist",
    "perc_mito_per_cell_hist",
    "highest_exprs_genes_dens_plt",
]


def plot_test(branch: CellBranch, **kwargs):
    run_name = branch.current_process
    plot_path = branch[run_name].plot_map["plot_test"]
    plt.plot([0, 1, 2, 3, 4], [0, 3, 5, 9, 11], **kwargs)
    plt.savefig(plot_path)


def qc_normalize(branch: "CellBranch", process: str, **kwargs):
    # allow parametrization of what QC plots to include
    _qc_plots = qc_plots if qc_plots in kwargs else DEFAULT_QC_PLOTS_NORMALIZE
    path_map = branch[process].path_map
    initial_process = branch.current_process if branch.current_process != None else "root"
    branch.goto_process(process)

    for qc_plot in _qc_plots:
        if qc_plot not in ALL_QC_PLOTS_NORMALIZE:
            raise ValueError(f"{qc_plot} QC plot does not exist. Select a subset from {ALL_QC_PLOTS_NORMALIZE}")

        locals()[qc_plot](branch, path_map[qc_plot])  # call on plotting functions

    branch.goto_process(initial_process)

