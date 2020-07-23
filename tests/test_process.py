from copy import deepcopy

import pytest
import pandas as pd

from cellforest import CellForest, Counts
from tests.fixtures import *
import tests
from tests.test_init import build_root_fix
from tests.test_data_ops import test_subset_fix

# TODO: add output file checks


@pytest.fixture
def test_normalize_fix(root_path, build_root_fix, norm_spec):
    cf = CellForest(root_dir=root_path, spec=norm_spec)
    cf.process.normalize()
    return cf


def test_process_chain(root_path, build_root_fix, process_chain_spec):
    cf = CellForest(root_dir=root_path, spec=process_chain_spec)
    cf.process.normalize()
    cf.process.test_process()
    return cf


def test_logging(test_normalize_fix):
    # TODO: QUEUE
    pass


def test_process_aliasing(root_path_2, sample_paths, alias_spec):
    cf = CellForest.from_input_dirs(root_path_2, sample_paths, spec=alias_spec, mode="rna")
    cf.process.process_1()
    cf.process.process_2()
    return cf


def test_normalize_cf_goto(test_subset_fix):
    cf = test_subset_fix
    rna = Counts.load(cf["normalize"].path_map["rna"])
    cf.goto_process("root")
    assert len(cf.meta) == 400
    assert len(cf.rna) == 400
    cf = cf.goto_process("normalize")
    assert len(cf.meta) == 11
    assert len(cf.rna) == 11
    assert cf.rna.shape == rna.shape
    assert len(cf.rna.features) == len(rna.features)


def test_reduce_on_existing_normalize():
    """
    Test that reduce works on an existing normalize directory and normalize
    isn't re-run
    """
    pass


def test_fresh_normalize_reduce():
    """
    Test that normalize and reduce work sequentially from scratch
    """
    pass
