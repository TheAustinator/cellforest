from cellforest import from_metadata
from tests.fixtures import *


def test_fork():
    pass


@pytest.fixture
def cb_1(norm_reduce_spec, merge_root_1, metadata):
    # add two test processes, which won't be run
    test_process_specs = [
        {"process": "test_process", "alias": "test_process_1",},
        {"process": "test_process", "alias": "test_process_2",},
    ]
    spec = norm_reduce_spec + test_process_specs
    cb = from_metadata(merge_root_1, metadata, spec=spec)
    cb.process.normalize()
    cb.process.reduce()
    return cb


def test_merge_empty(cb_1, merge_root_2):
    cb_1._merge(merge_root_2)
    # TODO: check dirs


def test_merge_partial_overlap():
    pass


def test_merge_no_overlap():
    pass


def test_merge_different_run_id():
    pass


def test_merge_multiple_subprocess_dirs():
    pass
