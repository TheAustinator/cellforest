from copy import deepcopy

import pandas as pd
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
def root_path_4(data_dir):
    return data_dir / "root_4"


@pytest.fixture
def root_path_5(data_dir):
    return data_dir / "root_5"


@pytest.fixture
def merge_root_1(data_dir):
    return data_dir / "merge/root_1"


@pytest.fixture
def merge_root_2(data_dir):
    return data_dir / "merge/root_2"


def root_path_example(data_dir):
    return data_dir / "example_usage" / "root"


@pytest.fixture
def sample_metadata_path(data_dir):
    return data_dir / "sample_metadata.tsv"


@pytest.fixture
def sample_metadata(sample_metadata_path):
    return pd.read_csv(sample_metadata_path, sep="\t")


@pytest.fixture
def counts_path(root_path):
    return root_path / "rna.pickle"


@pytest.fixture
def branch_spec_norm():
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
        },
    ]
    return spec


@pytest.fixture
def branch_spec_norm_reduce(branch_spec_norm):
    spec = deepcopy(branch_spec_norm)


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
        "_PROCESS_": "reduce",
        "_PARAMS_": {
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
def branch_spec_norm_reduce_cluster(branch_spec_norm_reduce):
    spec = deepcopy(branch_spec_norm_reduce)
    spec_run_cluster = {"_PROCESS_": "cluster", "_PARAMS_": {"num_pcs": 3, "res": 0.5, "eps": 0.1,}}
    spec.append(spec_run_cluster)
    return spec


@pytest.fixture
def process_chain_spec(branch_spec_norm):
    spec = deepcopy(branch_spec_norm)
    spec.append({"_PROCESS_": "test_process"})
    return spec


@pytest.fixture
def alias_spec():
    spec = [
        {"_PROCESS_": "test_process", "_ALIAS_": "process_1",},
        {"_PROCESS_": "test_process", "_ALIAS_": "process_2",},
    ]
    return spec


@pytest.fixture
def processes_of_norm_reduce_spec(norm_reduce_spec):
    avail_processes = []
    for process_spec in norm_reduce_spec:
        avail_processes.append(process_spec["process"])

    return avail_processes
