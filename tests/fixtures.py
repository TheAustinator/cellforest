from copy import deepcopy

import pandas as pd
from random import choice
from pathlib import Path
import pytest

from tests.utils.get_test_data import get_test_data


@pytest.fixture
def data_dir():
    path = Path(__file__).parent / "data"
    check_paths = [path / subdir for subdir in ["v2", "v3", "v3_gz"]]
    if not all(list(map(lambda p: p.exists(), check_paths))):
        get_test_data()
    return path


@pytest.fixture
def sample_1(data_dir):
    return data_dir / "v3/sample_1"


@pytest.fixture
def sample_2(data_dir):
    return data_dir / "v3/sample_2"


@pytest.fixture
def sample_1_v2(data_dir):
    return data_dir / "v2/sample_1"


@pytest.fixture
def sample_1_gz(data_dir):
    return data_dir / "v3_gz/sample_1"


@pytest.fixture
def sample_paths(sample_1, sample_2):
    return [sample_1, sample_2]


@pytest.fixture
def root_path(data_dir):
    return data_dir / "root_1"


@pytest.fixture
def root_path_2(data_dir):
    return data_dir / "root_2"


@pytest.fixture
def root_path_3(data_dir):
    return data_dir / "root_3"


@pytest.fixture
def root_path_example(data_dir):
    return data_dir / "example_usage" / "root"


@pytest.fixture
def metadata_path(data_dir):
    return data_dir / "sample_metadata.tsv"


@pytest.fixture
def metadata(metadata_path):
    return pd.read_csv(metadata_path, sep="\t")


@pytest.fixture
def counts_path(root_path):
    return root_path / "rna.pickle"


@pytest.fixture
def norm_spec():
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
        }
    ]
    return spec


@pytest.fixture
def norm_sctransform_spec():
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
        }
    ]
    return spec


@pytest.fixture
def norm_reduce_spec(norm_spec):
    spec = deepcopy(norm_spec)
    reduce_run_spec = {
        "process": "reduce",
        "params": {
            "pca_npcs": 3,
            "umap_n_neighbors": 3,
            "umap_min_dist": 0.1,
            "umap_n_components": 2,
            "umap_metric": "euclidean",
        },
    }
    spec.append(reduce_run_spec)
    return spec


@pytest.fixture
def process_chain_spec(norm_spec):
    spec = deepcopy(norm_spec)
    spec.append({"process": "test_process"})
    return spec


@pytest.fixture
def alias_spec():
    spec = [{"process": "test_process", "alias": "process_1",}, {"process": "test_process", "alias": "process_2",}]
    return spec


@pytest.fixture
def random_process():
    AVAILABLE_PROCESSES = ["root", "normalize", "reduce"]
    return choice(AVAILABLE_PROCESSES)
