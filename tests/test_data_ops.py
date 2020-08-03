import cellforest
from tests.test_init import *


@pytest.fixture
def test_subset_fix(root_path_3, metadata, norm_spec):
    spec = deepcopy(norm_spec)
    spec[0]["_SUBSET_"] = {"sample": "sample_1"}
    branch = cellforest.from_metadata(root_path_3, metadata, spec=spec)
    branch.process.normalize()
    output_meta_path = branch["normalize"].path_map["meta"]
    output_meta = pd.read_csv(output_meta_path, sep="\t")
    assert (output_meta["sample"] == "sample_1").all()
    assert len(branch.meta) == len(branch.rna)
    return branch


def test_subset(test_subset_fix):
    pass


def test_subset_multiple(root_path_3, metadata, norm_spec):
    pass


def test_filter(build_root_fix):
    pass


def test_partition(build_root_fix):
    pass


def goto(test_subset_fix):
    branch = test_subset_fix
    branch.goto("root")
