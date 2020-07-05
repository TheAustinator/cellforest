from cellforest import Counts
from tests.fixtures import *


@pytest.fixture
def test_from_cellranger_fix(sample_1):
    rna = Counts.from_cellranger(sample_1)
    return rna


@pytest.fixture
def test_save_fix(test_from_cellranger_fix, counts_path):
    test_from_cellranger_fix.save(counts_path)
    return counts_path


def test_from_cellranger_v2(sample_1_v2):
    rna = Counts.from_cellranger(sample_1_v2)
    return rna


def test_from_cellranger_gz(sample_1_gz):
    rna = Counts.from_cellranger(sample_1_gz)
    return rna


def test_load(test_save_fix):
    rna = Counts.load(test_save_fix)
    return rna


def test_concatenate(test_from_cellranger_fix):
    rna = test_from_cellranger_fix[:50, :50]
    assert rna.append(rna).shape[0] == 100
    assert rna.append(rna, axis=1).shape[1] == 100
    assert rna.append([rna, rna]).shape[0] == 150


def test_slice(test_from_cellranger_fix):
    rna = test_from_cellranger_fix
    assert rna.shape == (rna.cell_ids.shape[0], rna.genes.shape[0])
    assert rna["AAACATACAACCAC-1"].shape[0] == 1
    assert rna["AAACATACAACCAC-1", :].shape[0] == 1
    assert rna[["AAACATACAACCAC-1", "AAACATTGATCAGC-1"], :].shape[0] == 2
    assert rna[pd.Series(["AAACATACAACCAC-1", "AAACATTGATCAGC-1"]), :].shape[0] == 2
    assert rna[:, "ENSG00000243485"].shape[1] == 1
    assert rna[:, ["ENSG00000243485", "ENSG00000186092"]].shape[1] == 2
    assert rna[:, pd.Series(["ENSG00000243485", "ENSG00000186092"])].shape[1] == 2
    assert rna[:, "FAM138A"].shape[1] == 1
    assert rna[:, ["OR4F5", "FAM138A"]].shape[1] == 2
    assert rna[:, pd.Series(["OR4F5", "FAM138A"])].shape[1] == 2
    assert rna[1, :].shape[0] == 1
    assert rna[:, 1].shape[1] == 1
    assert rna[5:10, 5:10].shape == (5, 5)


def test_from_cellranger(test_from_cellranger_fix):
    pass


def test_save(test_save_fix):
    pass
