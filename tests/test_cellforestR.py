import pytest
from re import escape

from cellforest.utils.r.run_r_script import run_r_script
from tests.fixtures import *

# TODO: Add this function to test_quick.py
def test_cellforest_load_last_process(
    root_path_example, norm_reduce_spec, processes_of_norm_reduce_spec, test_norm_reduce
):
    path_to_cellforestR_script = Path(__file__).parent / "r" / "cellforest_load.R"
    # TODO: Add support for passing spec through Rscript
    run_r_script(path_to_cellforestR_script, [root_path_example, "placeholder", processes_of_norm_reduce_spec[-1]])


# TODO: Add this function to test_all.py
def test_cellforest_load_all_processes(
    root_path_example, norm_reduce_spec, processes_of_norm_reduce_spec, test_norm_reduce
):
    path_to_cellforestR_script = Path(__file__).parent / "r" / "cellforest_load.R"
    # TODO: Add support for passing spec through Rscript
    for process in processes_of_norm_reduce_spec:
        run_r_script(path_to_cellforestR_script, [root_path_example, "placeholder", process])
