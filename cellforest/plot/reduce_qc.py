import matplotlib
import matplotlib.pyplot as plt

from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_process_r_script
from cellforest.plot.qc_plot_wrappers import qc_plot_py, qc_plot_r


@qc_plot_r
def plot_pca_elbow_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_process_r_script(branch, r_script, args, branch.current_process)


@qc_plot_r
def plot_pca_loadings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_process_r_script(branch, r_script, args, branch.current_process)


@qc_plot_r
def plot_pca_embeddings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_process_r_script(branch, r_script, args, branch.current_process)


@qc_plot_r
def plot_umap_embeddings_scat(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_process_r_script(branch, r_script, args, branch.current_process)
