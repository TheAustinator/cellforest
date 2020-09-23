from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_r_script_logged
from cellforest.plot.qc_plot_wrappers import qc_plot_r


@qc_plot_r
def plot_pca_elbow_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_elbow_curv")


@qc_plot_r
def plot_pca_loadings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_loadings_scat")


@qc_plot_r
def plot_pca_embeddings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_pca_embeddings_scat")


@qc_plot_r
def plot_umap_embeddings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_umap_embeddings_scat")


@qc_plot_r
def plot_cell_cycle_scoring_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, "plot_cell_cycle_scoring_scat")
