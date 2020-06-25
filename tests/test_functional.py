import pytest
import pandas as pd

from cellforest import CellForest


def test_init(root_dir):
    spec = {
        "combine": {},
        "normalize": {
            "min_genes": 5,
            "max_genes": 5000,
            "min_cells": 5,
            "prec_mito_cutoff": 20,
            "method": "seurat_default",
        },
    }
    cf = CellForest(root_dir=root_dir, spec_dict=spec)
    assert cf.rna["TTTGTCAGTTTAGCTG-1"].shape[0] == 1
    assert cf.rna["TTTGTCAGTTTAGCTG-1", :].shape[0] == 1
    assert cf.rna[["TTTGTCAGTTTAGCTG-1", "TTTGTCAGTTCAACCA-1"], :].shape[0] == 2
    assert cf.rna[pd.Series(["TTTGTCAGTTTAGCTG-1", "TTTGTCAGTTCAACCA-1"]), :].shape[0] == 2
    assert cf.rna[:, "ENSG00000268674"].shape[1] == 1
    assert cf.rna[:, ["ENSG00000268674", "ENSG00000238009"]].shape[1] == 2
    assert cf.rna[:, pd.Series(["ENSG00000268674", "ENSG00000238009"])].shape[1] == 2
    assert cf.rna[:, "FAM138A"].shape[1] == 1
    assert cf.rna[:, ["RP11-34P13.3", "FAM138A"]].shape[1] == 2
    assert cf.rna[:, pd.Series(["RP11-34P13.3", "FAM138A"])].shape[1] == 2
    assert cf.rna[1, :].shape[0] == 1
    assert cf.rna[:, 1].shape[1] == 1
    assert cf.rna[5:10, 5:10].shape == (5, 5)


if __name__ == "__main__":
    test_init("/Users/austinmckay/data")
