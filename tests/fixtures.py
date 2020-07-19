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
def metadata_path(data_dir):
    return data_dir / "sample_metadata.tsv"


@pytest.fixture
def metadata(metadata_path):
    return pd.read_csv(metadata_path, sep="\t")


@pytest.fixture
def counts_path(root_path):
    return root_path / "rna.pickle"
