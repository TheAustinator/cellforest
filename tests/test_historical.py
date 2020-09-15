import cellforest as cf
from tests.fixtures import *


def test_spec_change(data_dir, sample_metadata, root_path_4):
    root = root_path_4
    spec = [
        {
            "_PROCESS_": "normalize",
            "_PARAMS_": {
                "min_genes": 5,
                "max_genes": 5000,
                "min_cells": 5,
                "nfeatures": 30,
                "perc_mito_cutoff": 20,
                "method": "seurat_default",
            },
            "_SUBSET_": {"sample_id": "sample_1"},
        }
    ]
    branch = cf.from_sample_metadata(root, sample_metadata, branch_spec=spec)
    branch.process.normalize()
    branch.rna  # works fine during first go
    spec = [
        {
            "_PROCESS_": "normalize",
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
    branch = cf.load(root, branch_spec=spec)
    branch.process.normalize()  # breaks
    branch.rna  # breaks
