import pandas as pd
from pathlib import Path
import pytest

from cellforest import Counts
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
def metadata(data_dir):
    path = data_dir / "sample_metadata.tsv"
    base_path = Path(__file__).parent / "data/v3"
    cellranger_paths = [base_path / "sample_1", base_path / "sample_2"]
    if not path.exists():
        df = pd.DataFrame({"name": ["sample_1", "sample_2"], "path": cellranger_paths})
        df.to_csv(path, sep="\t")
    else:
        df = pd.read_csv(path, sep="\t")
    return df


@pytest.fixture
def counts_path(root_path):
    return root_path / "counts.pickle"


@pytest.fixture
def test_from_cellranger(sample_1):
    rna = Counts.from_cellranger(sample_1)
    return rna


@pytest.fixture
def test_save(test_from_cellranger, counts_path):
    test_from_cellranger.save(counts_path)
    return counts_path
