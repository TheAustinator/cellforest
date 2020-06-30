from cellforest import CellForest
from tests.fixtures import *


def test_from_input_dirs(root_path, sample_paths):
    cf = CellForest.from_input_dirs(root_path, sample_paths, mode="rna")
    _ = cf.meta
    _ = cf.rna
    return cf


def test_from_metdata(root_path, metadata):
    cf = CellForest.from_metadata(root_path, metadata)
    _ = cf.meta
    _ = cf.rna
    assert len(cf.meta.columns) > 0
    return cf
