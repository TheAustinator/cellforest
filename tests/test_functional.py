import pytest
from cellforest import CellForest


def test_init(root_dir):
    spec = {}
    cf = CellForest(root_dir=root_dir, spec_dict=spec)
    rna = cf.rna
    print(rna)


if __name__ == "__main__":
    test_init("/Users/austinmckay/data")
