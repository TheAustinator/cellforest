import pytest
from re import escape

from cellforest.utils.r.run_r_script import run_r_script
from tests.fixtures import *


def test_cellforestR(root_path_example, norm_spec, random_process):
    path_to_cellforestR_script = Path(__file__).parent / "r" / "cellforest_load.R"
    # TO-DO: Add support for passing spec through Rscript
    run_r_script(path_to_cellforestR_script, [root_path_example, "placeholder", random_process])
