import pytest
import pandas as pd

from cellforest import CellForest, Counts
from tests.fixtures import *
import tests
from tests.test_init import build_root_fix


@pytest.fixture
def test_normalize_fix(root_path, build_root_fix):
    spec = {
        "normalize": {
            "min_genes": 5,
            "max_genes": 5000,
            "min_cells": 5,
            "nfeatures": 30,
            "perc_mito_cutoff": 20,
            "method": "seurat_default",
        },
    }
    cf = CellForest(root_dir=root_path, spec_dict=spec)
    cf.process.normalize()
    return cf


def test_logging(test_normalize_fix):
    # TODO: QUEUE
    pass


def test_normalize_cf_at(test_normalize_fix):
    """Functionality not yet implemented"""
    return

    cf = test_normalize_fix
    rna = Counts.load(cf["normalize"].path_map["rna"])
    cf = cf.at("normalize")
    assert cf.rna.shape == rna.shape
    assert len(cf.meta) == len(rna)
    assert len(cf.rna.features) == len(rna.features)
