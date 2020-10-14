from cellforest import CellBranch
from cellforest.utils.r.run_r_script import run_r_script_logged
from cellforest.plot.qc_plot_wrappers import qc_plot_py, qc_plot_r


@qc_plot_r
def plot_umis_per_barcode_rank_curv(branch: "CellBranch", r_script: str, args: list, **kwargs):
    run_r_script_logged(branch, r_script, args, branch.current_process)

