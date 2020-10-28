from copy import deepcopy
from typing import Type


class DefaultsMeta(type):
    SPEC_BASE = [
        {
            "_PROCESS_": "normalize",
            "_PARAMS_": {
                "min_genes": 0,
                "max_genes": 37000,
                "min_cells": 3,
                "nfeatures": 3000,
                "perc_mito_cutoff": 100,
                "method": "seurat_default",
            },
        },
        {
            "_PROCESS_": "reduce",
            "_PARAMS_": {
                "pca_npcs": 30,
                "umap_n_neighbors": 30,
                "umap_min_dist": 0.3,
                "umap_n_components": 2,
                "umap_metric": "correlation",
            },
        },
        {"_PROCESS_": "cluster", "_PARAMS_": {"num_pcs": 3, "res": 0.5, "eps": 0.1,}},
    ]

    SPEC_MARKERS = deepcopy(SPEC_BASE) + [
        {"_PROCESS_": "markers", "_PARAMS_": {"logfc_thresh": 0.25, "test": "wilcox"}},
    ]

    SPEC_DIFFEXP = deepcopy(SPEC_MARKERS) + [
        {"_PROCESS_": "diffexp", "_PARAMS_": {"logfc_thresh": 0.25, "test": "wilcox"}, "_PARTITION_": {"condition"}},
    ]

    @property
    def spec_markers(cls):
        return deepcopy(cls.SPEC_MARKERS)

    @property
    def spec_diffexp(cls):
        return deepcopy(cls.SPEC_DIFFEXP)


class defaults(metaclass=DefaultsMeta):
    pass


def build_const_container(container_name: str, const_dict: dict):
    container_meta = type(f"{container_name}Meta", {"__init__": _container_meta_init})
    for k, v in const_dict.items():
        setattr(container_meta, k, v)
    # container =


class ConstContainer:
    def __init__(self, container_name: str, const_dict: dict):
        for k, v in const_dict:
            setattr(self, k, v)

    @classmethod
    def build(cls, container_name: str, const_dict: dict) -> Type:
        raise NotImplementedError()
