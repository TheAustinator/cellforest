import pytest

import cellforest as cf
from tests.fixtures import *


def test_from_input_dirs_single(root_path, sample_1):
    branch = cf.from_input_dirs(root_path, sample_1, mode="rna")
    _ = branch.meta
    _ = branch.rna
    return sample_1


def test_from_input_dirs(test_from_input_dirs_fix):
    pass


def test_from_meta(build_root_fix):
    pass
