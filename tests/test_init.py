import pytest

from cellforest import CellBranch
from tests.fixtures import *


@pytest.fixture
def test_from_input_dirs_fix(root_path, sample_paths):
    cf = CellBranch.from_input_dirs(root_path, sample_paths, mode="rna")
    _ = cf.meta
    _ = cf.rna
    return sample_paths


@pytest.fixture
def build_root_fix(root_path, metadata):
    cf = CellBranch.from_metadata(root_path, metadata)
    _ = cf.meta
    _ = cf.rna
    assert len(cf.meta.columns) > 0
    return cf


def test_from_input_dirs_single(root_path, sample_1):
    cf = CellBranch.from_input_dirs(root_path, sample_1, mode="rna")
    _ = cf.meta
    _ = cf.rna
    return sample_1


def test_from_input_dirs(test_from_input_dirs_fix):
    pass


def test_from_meta(build_root_fix):
    pass
