import pytest


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
                "pca_npcs": {"_SWEEP_": {"min": 4, "max": 5, "base": 2}},
                "umap_n_neighbors": {"_SWEEP_": {"min": 2, "max": 3, "step": 1}},
                "umap_min_dist": 0.1,
                "umap_n_components": 2,
                "umap_metric": {"_SWEEP_": ["cosine", "euclidean"]},
            },
        },
    ]
    return spec


def test_datatree():
    tree =