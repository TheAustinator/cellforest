import pytest

import cellforest as cf
from tests.fixtures import *


@pytest.fixture
def tree_spec():
    # for a given process, either matched arrays for everything or  that sweep everything, or
    spec = [
        {
            "_PROCESS_": "normalize",
            "_PARAMS_": [
                {
                    "min_genes": 5,
                    "max_genes": {"_SWEEP_": {"min": 2000, "max": 3000, "step": 1000}},
                    "min_cells": 5,
                    "nfeatures": {"_SWEEP_": {"min": 2, "max": 3, "base": 2}},
                    "perc_mito_cutoff": 20,
                    "method": "seurat_default",
                },
                {"min_genes": 5, "max_genes": 5000, "min_cells": 5, "perc_mito_cutoff": 20, "method": "sctransform"},
            ],
            "_SUBSET_": {"sample": {"_SWEEP_": ["sample_1", "sample_2"]}},
        },
        {
            "_PROCESS_": "reduce",
            "_PARAMS_": {
                "pca_npcs": 16,
                "umap_n_neighbors": 3,
                "umap_min_dist": 0.1,
                "umap_n_components": 2,
                "umap_metric": "euclidean",
            },
        },
    ]
    return spec


@pytest.fixture
def tree_spec_reduce_update():
    spec = {
        "_PROCESS_": "reduce",
        "_PARAMS_": {
            "pca_npcs": {"_SWEEP_": {"min": 4, "max": 5, "base": 2}},
            "umap_n_neighbors": {"_SWEEP_": {"min": 2, "max": 3, "step": 1}},
            "umap_min_dist": 0.1,
            "umap_n_components": 2,
            "umap_metric": {"_SWEEP_": ["cosine", "euclidean"]},
        },
    }
    return spec


@pytest.fixture
def tree_spec_norm_update():
    spec = {
        "_PROCESS_": "normalize",
        "_PARAMS_": {
            "min_genes": 5,
            "max_genes": 2000,
            "min_cells": 5,
            "nfeatures": 4,
            "perc_mito_cutoff": 20,
            "method": "seurat_default",
        },
        "_SUBSET_": {"sample": "sample_1"},
    }
    return spec


def test_datatree(root_path_5, sample_metadata, tree_spec, tree_spec_reduce_update, tree_spec_norm_update):
    tree = cf.from_sample_metadata(root=root_path_5, sample_metadata=sample_metadata, tree_spec=tree_spec)
    assert tree.n_branches == 10
    tree.process.normalize()
    tree.process.reduce()
    tree.update_process_spec("reduce", tree_spec_reduce_update)
    assert tree.n_branches == 80
    tree.update_process_spec("normalize", tree_spec_norm_update)
    assert tree.n_branches == 8
    tree.process.reduce()
    spec_str = str(tree.tree_spec.branch_specs[0])
    reduce_sweep = tree._branch_cache[spec_str]["normalize"].path / "reduce"
    assert len(list(reduce_sweep.iterdir())) >= 8
    return tree
