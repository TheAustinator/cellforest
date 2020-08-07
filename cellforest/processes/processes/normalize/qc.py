from pathlib import Path
import matplotlib.pyplot as plt

from dataforest.hooks import dataprocess

from cellforest.utils.r.run_r_script import run_process_r_script

ALL_QC_PLOTS_NORMALIZE = [  # reference of all supported QC plots
    "genes_per_cell_hist",
    "umis_per_cell_hist",
    "umis_vs_genes_hist",
    "perc_mito_per_cell_hist",
    "highest_expression_genes_dens_plt",
]

DEFAULT_QC_PLOTS_NORMALIZE = [  # default QC plots
    "genes_per_cell_hist",
    "umis_per_cell_hist",
    "umis_vs_genes_hist",
    "perc_mito_per_cell_hist",
    "highest_expression_genes_dens_plt",
]


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


def genes_per_cell_hist(branch: "CellBranch", save_path: str):
    branch.rna.hist("nonzero", axis=0)
    plt.savefig(save_path, dpi=330)


def umis_per_cell_hist(branch: "CellBranch", save_path: str):
    branch.rna.hist("sum", axis=0)
    plt.savefig(save_path, dpi=330)


def umis_vs_genes_hist(branch: "CellBranch", save_path: str):
    branch.rna.scatter(agg_x="nonzero", agg_y="sum", axis=0)


def perc_mito_per_cell_hist(branch: "CellBranch", save_path: str):
    branch.rna.hist()


def highest_expression_genes_dens_plt(branch: "CellBranch", save_path: str):
    pass
    # R script
