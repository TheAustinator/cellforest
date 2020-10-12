from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_r_script_logged
from cellforest.plot.qc_plot_wrappers import qc_plot_r


@qc_plot_r
def plot_marker_genes_per_cluster_bar(branch: "CellBranch", r_script: str, args: list, **kwargs):
    args.append(branch.current_path)
    run_r_script_logged(branch, r_script, args, "plot_marker_genes_per_cluster_bar")
