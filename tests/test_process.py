import pytest

from cellforest import CellBranch, Counts
from tests.fixtures import *
import tests
from tests.test_init import build_root_fix
from tests.test_data_ops import test_subset_fix

# TODO: add output file checks


@pytest.fixture
def test_norm_fix(root_path, build_root_fix, norm_spec):
    cf = CellBranch(root_dir=root_path, spec=norm_spec)
    cf.process.normalize()
    return cf


def test_norm_reduce(root_path, build_root_fix, norm_reduce_spec, test_norm_fix):
    cf = CellBranch(root_dir=root_path, spec=norm_reduce_spec)
    cf.process.reduce()
    return cf


def test_process_chain(root_path, build_root_fix, process_chain_spec):
    cf = CellBranch(root_dir=root_path, spec=process_chain_spec)
    cf.process.normalize()
    cf.process.test_process()
    return cf


def test_logging(test_norm_fix):
    # TODO: QUEUE
    pass


def test_process_aliasing(root_path_2, sample_paths, alias_spec):
    cf = CellBranch.from_input_dirs(root_path_2, sample_paths, spec=alias_spec, mode="rna")
    cf.process.process_1()
    cf.process.process_2()
    return cf


def test_normalize_cf_goto(test_subset_fix):
    cf = test_subset_fix
    rna = Counts.load(cf["normalize"].path_map["rna"])
    cf.goto_process("root")
    assert len(cf.meta) == 600
    assert len(cf.rna) == 600
    cf = cf.goto_process("normalize")
    assert len(cf.meta) == 59
    assert len(cf.rna) == 59
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
