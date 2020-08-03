import cellforest
from tests.fixtures import *


def test_spec_change(data_dir, metadata, root_path_4):
    root = root_path_4
    spec = [
        {
            "process": "normalize",
            "_PARAMS_": {
                "min_genes": 5,
                "max_genes": 5000,
                "min_cells": 5,
                "nfeatures": 30,
                "perc_mito_cutoff": 20,
                "method": "seurat_default",
            },
            "_SUBSET_": {"sample": "sample_1"},
        }
    ]
    branch = cellforest.from_metadata(root, metadata, spec=spec)
    branch.process.normalize()
    branch.rna  # works fine during first go
    spec = [
        {
            "process": "normalize",
            "_PARAMS_": {
                "min_genes": 2,
                "max_genes": 5000,
                "min_cells": 5,
                "nfeatures": 30,
                "perc_mito_cutoff": 20,
                "method": "seurat_default",
            },
            "_SUBSET_": {"sample": "sample_1"},
        }
    ]
    branch = cellforest.load(root, spec=spec)
    branch.process.normalize()  # breaks
    branch.rna  # breaks
