import pytest


@pytest.fixture
def tree_spec():
    # for a given process, either matched arrays for everything or  that sweep everything, or
    tree_spec = [{"process": "normalize", "_PARAMS_": [{}, {}], "_SUBSET_": [{""}]}, {}]
