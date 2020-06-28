import pandas as pd
from pathlib import Path
import pytest

from tests.utils.get_test_data import get_test_data


@pytest.fixture
def data_dir():
    return Path(__file__).parent / "data"


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
def metadata_path(data_dir):
    path = data_dir / "sample_metadata.tsv"
    if not path.exists:
        df = pd.DataFrame({"name": ["sample_1", "sample_2"], "path": sample_paths})
        df.to_csv(path, sep="\t")
    return metadata_path


@pytest.fixture
def counts_path(data_dir):
    return data_dir / "counts.pickle"
