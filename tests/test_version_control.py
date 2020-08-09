import cellforest as cf
from tests.fixtures import *


def test_fork():
    pass


@pytest.fixture
def cb_1(branch_spec_norm_reduce, merge_root_1, sample_metadata):
    # add two test processes, which won't be run
    test_process_specs = [
        {"_PROCESS_": "test_process", "_ALIAS_": "test_process_1",},
        {"_PROCESS_": "test_process", "_ALIAS_": "test_process_2",},
    ]
    spec = branch_spec_norm_reduce + test_process_specs
    cb = cf.from_sample_metadata(merge_root_1, sample_metadata, branch_spec=spec)
    cb.process.normalize()
    cb.process.reduce()
    return cb


def test_merge_empty(cb_1, merge_root_2):
    # TODO: need to build merge_root_2 first
    # cb_1._merge(merge_root_2)
    # TODO: check dirs
    pass


def test_merge_partial_overlap():
    pass


def test_merge_no_overlap():
    pass


def test_merge_different_run_id():
    pass


def test_merge_multiple_subprocess_dirs():
    pass
