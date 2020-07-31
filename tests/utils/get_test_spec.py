def get_spec():
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


def get_spec_1():
    spec = [
        {
            "process": "normalize",
            "params": {
                "min_genes": 3,
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


def get_spec_2():
    spec = [
        {
            "process": "normalize",
            "params": {
                "min_genes": 5,
                "max_genes": 5000,
                "min_cells": 5,
                "perc_mito_cutoff": 20,
                "method": "sctransform",
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


def get_spec_3():
    spec = [
        {
            "process": "normalize",
            "params": {
                "min_genes": 4,
                "max_genes": 5001,
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
