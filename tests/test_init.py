import pytest

import cellforest as cf
from tests.fixtures import *


@pytest.fixture
def test_from_input_dirs_fix(root_path, sample_paths):
    branch = cf.from_input_dirs(root_path, sample_paths, mode="rna")
    _ = branch.meta
    _ = branch.rna
    return sample_paths


@pytest.fixture
def build_root_fix(root_path, sample_metadata):
    branch = cf.from_sample_metadata(root_path, sample_metadata)
    _ = branch.meta
    _ = branch.rna
    assert len(branch.meta.columns) > 0
    return branch


def test_from_input_dirs_single(root_path, sample_1):
    branch = cf.from_input_dirs(root_path, sample_1, mode="rna")
    _ = branch.meta
    _ = branch.rna
    return sample_1


def test_from_input_dirs(test_from_input_dirs_fix):
    pass


def test_from_meta(build_root_fix):
    pass
