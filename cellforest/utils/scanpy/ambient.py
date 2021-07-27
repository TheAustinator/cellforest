import os
from pathlib import Path
from typing import AnyStr, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from cellforest.utils import r
from cellforest.utils.shell.shell_command import process_shell_command

_R_UTILS_DIR = Path(r.__file__).parent
_R_RUN_AMBIENT_RNA = str(_R_UTILS_DIR / "run_ambient_rna.R")


def est_ambient_rna(input_dir_10x: AnyStr, output_dir: AnyStr, input_path_clusters: Optional[AnyStr] = None):
    output_dir = Path(output_dir)
    proc_dir = output_dir / "process_files"
    os.makedirs(proc_dir, exist_ok=True)
    input_path_clusters = input_path_clusters if input_path_clusters else ""
    output_path_dx_est = str(proc_dir / "decontx_est.csv")
    output_path_dx_h5 = str(output_dir / "decontx_counts.h5")
    output_path_sx_est = str(proc_dir / "soupx_est.csv")
    output_path_sx_prof = str(output_dir / "soupx_prof.csv")
    output_path_sx_h5 = str(output_dir / "soupx_counts.h5")
    arg_list = [
        input_dir_10x,
        input_path_clusters,
        output_path_dx_est,
        output_path_dx_h5,
        output_path_sx_est,
        output_path_sx_prof,
        output_path_sx_h5,
    ]
    command_string = f"Rscript {_R_RUN_AMBIENT_RNA} {' '.join(map(str, arg_list))}"
    process_shell_command(command_string=command_string, logs_dir=proc_dir, logfile_prefix="est_ambient_rna")

    dx_est = pd.read_csv(output_path_dx_est, index_col=0).set_index("Barcode").drop(columns="Sample")
    dx_est.columns = dx_est.columns.str.lower()
    soupx_ran = all([os.path.exists(p) for p in [output_path_sx_prof, output_path_sx_est]])
    if soupx_ran:
        sx_prof = pd.read_csv(output_path_sx_prof, index_col=0)
        sx_est = pd.read_csv(output_path_sx_est, index_col=0).rename(columns={"rho": "soupx_contamination"})
        sx_prof.to_csv(output_path_sx_prof)
        sx_est.to_csv(output_path_sx_est)
        df_est = dx_est.merge(sx_est, left_index=True, right_index=True)
    else:
        df_est = dx_est
        df_est["soupx_contamination"] = np.nan
    dx_est.to_csv(output_path_dx_est)
    df_est.to_csv(output_dir / "cell_contamination.csv")
    df_est[["decontx_contamination", "soupx_contamination"]].mean()
    est_mean = df_est[["decontx_contamination", "soupx_contamination"]].mean()
    est_90th = df_est[["decontx_contamination", "soupx_contamination"]].apply(lambda arr: np.percentile(arr, 90))
    est_summary = pd.DataFrame({"mean": est_mean, "90th": est_90th})
    est_summary.to_csv(output_dir / "sample_contamination.csv")
    # plot
    bins = np.arange(0, 1.025, 0.025)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 3))
    ax1.hist(df_est["decontx_contamination"], label="decontx", bins=bins, log=True, alpha=0.25)
    if soupx_ran:
        ax1.hist(df_est["soupx_contamination"], label="soupx", bins=bins, log=True, alpha=0.25)
        ax2.hist(sx_prof["est"], bins=bins, log=True, alpha=0.5)
    ax1.set_title("cell contamination")
    ax2.set_title("gene contamination")
    ax1.set_xlabel("contamination")
    ax2.set_xlabel("contamination")
    ax1.set_ylabel("cells")
    ax2.set_ylabel("genes")
    ax1.legend()
    fig.savefig(output_dir / "contamination_plot.png")
