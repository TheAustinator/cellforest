from pathlib import Path
import matplotlib.pyplot as plt

from dataforest.hooks import dataprocess

# TODO: what to do about core/utility methods? core module? move to utils?
from cellforest.utils.r.run_r_script import run_process_r_script
from cellforest.processes import R_FUNCTIONS_FILEPATH


def qc_normalize(branch: "CellBranch", **kwargs):
    pass

