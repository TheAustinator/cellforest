import pytest

import cellforest as cf
from cellforest import CellBranch, Counts
from tests.fixtures import *
import tests
from tests.test_init import build_root_fix
from tests.test_data_ops import test_subset_fix

# TODO: add output file checks


@pytest.fixture
def test_norm_fix(root_path, build_root_fix, branch_spec_norm):
    branch = CellBranch(root=root_path, branch_spec=branch_spec_norm)
    branch.process.normalize()
    return branch


def test_norm_reduce(root_path, build_root_fix, branch_spec_norm_reduce, test_norm_fix):
    branch = CellBranch(root=root_path, branch_spec=branch_spec_norm_reduce)
    branch.process.reduce()
    return branch


def test_process_chain(root_path, build_root_fix, process_chain_spec):
    branch = CellBranch(root=root_path, branch_spec=process_chain_spec)
    branch.process.normalize()
    branch.process.test_process()
    return branch


def test_logging(test_norm_fix):
    # TODO: QUEUE
    pass


def test_process_aliasing_and_plotting(root_path_2, sample_paths, alias_spec):
    branch = cf.from_input_dirs(root_path_2, sample_paths, branch_spec=alias_spec, mode="rna")
    branch.process.process_1()
    branch.process.process_2()
    assert branch["process_2"].plot_map["plot_test"].exists()
    return branch


def test_normalize_branch_goto(test_subset_fix):
    branch = test_subset_fix
    rna = Counts.load(branch["normalize"].path_map["rna"])
    branch.goto_process("root")
    assert len(branch.meta) == 600
    assert len(branch.rna) == 600
    branch = branch.goto_process("normalize")
    assert len(branch.meta) == 59
    assert len(branch.rna) == 59
    assert branch.rna.shape == rna.shape
    assert len(branch.rna.features) == len(rna.features)


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
