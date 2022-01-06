import pytest

import cellforest as cf
from tests.fixtures import *
from tests.utils.check_files import check_root_files


def test_from_input_dirs_single(root_path, sample_1):
    branch = cf.from_input_dirs(root_path, sample_1, mode="rna", parallel=False)
    assert "AAACATACAACCAC-1" in branch.meta.index
    assert branch.rna.shape[0] > 1 and branch.rna.shape[1] > 1
    check_root_files(branch)
    return sample_1


def test_from_input_dirs(test_from_input_dirs_fix):
    pass


def test_from_meta(build_root):
    pass
