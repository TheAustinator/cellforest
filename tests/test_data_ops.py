import cellforest as cf
from tests.test_init import *


@pytest.fixture
def test_subset_fix(root_path_3, sample_metadata, branch_spec_norm):
    branch_spec = deepcopy(branch_spec_norm)
    branch_spec[0]["_SUBSET_"] = {"entity_id": "sample_1"}
    branch = cf.from_sample_metadata(root_path_3, sample_metadata, branch_spec=branch_spec)
    branch.process.normalize()
    output_meta_path = branch["normalize"].path_map["meta"]
    output_meta = pd.read_csv(output_meta_path, sep="\t")
    assert (output_meta["entity_id"] == "sample_1").all()
    assert len(branch.meta) == len(branch.rna)
    return branch


@useconfig("no_plot_config")
def test_subset(test_subset_fix):
    pass


def test_subset_multiple(root_path_3, sample_metadata, branch_spec_norm):
    pass


# TODO: filter and partition
# def test_filter(build_root):
#     pass


# def test_partition(build_root):
#     pass


def goto(test_subset_fix):
    branch = test_subset_fix
    branch.goto("root")
