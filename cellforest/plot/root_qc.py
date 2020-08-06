import matplotlib.pyplot as plt

from cellforest import CellBranch


def plot_test(branch: CellBranch, **kwargs):
    run_name = branch.current_process
    plot_path = branch[run_name].plot_map["plot_test"]
    plt.plot([0, 1, 2, 3, 4], [0, 3, 5, 9, 11], **kwargs)
    plt.savefig(plot_path)
