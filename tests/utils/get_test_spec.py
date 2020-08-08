def get_example_spec():
    spec = [
        {
            "process": "normalize",
            "params": {
                "min_genes": 5,
                "max_genes": 5000,
                "min_cells": 5,
                "nfeatures": 30,
                "perc_mito_cutoff": 20,
                "method": "seurat_default",
            },
            "subset": {"sample": "sample_1"},
        },
        {
            "process": "reduce",
            "params": {
                "pca_npcs": 3,
                "umap_n_neighbors": 3,
                "umap_min_dist": 0.1,
                "umap_n_components": 2,
                "umap_metric": "euclidean",
            },
        },
    ]
    return spec
