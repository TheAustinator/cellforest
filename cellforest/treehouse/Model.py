import numpy as np
import pandas as pd

from dataforest.treehouse.model import DataFrameModel


class CellModel(DataFrameModel):
    def __init__(self, tree):
        self.tree = tree
        super().__init__(self.tree.meta)

    @property
    def branch(self):
        # TODO: make this dynamic via menu
        return self.tree.branches[0]

    @property
    def df(self):
        return self.branch.meta

    @property
    def rna(self):
        return self.branch.rna

    @staticmethod
    def get_expr_df(rna, cell_ids, density=True):
        rna = rna[cell_ids]
        rna_avg = rna.sum(axis=0)
        rna_avg = np.squeeze(np.asarray(rna_avg))
        if density:
            rna_avg /= rna_avg.sum() / 100
        genes = rna.genes
        df = pd.DataFrame({"gene": genes, "expr": rna_avg.flatten()})
        df.sort_values("expr", ascending=False, inplace=True)
        df.index = df["gene"]
        return df

    def get_crude_markers(self, cells_1, cells_2):
        expr_1 = self.get_expr_df(self.rna, cells_1)
        expr_2 = self.get_expr_df(self.rna, cells_2)
        df = expr_1.copy()
        rna_root = self.branch.copy().goto_process("root").rna
        df["selected"] = self.get_expr_df(rna_root, cells_1, density=False)["expr"]
        df["rest"] = self.get_expr_df(rna_root, cells_2, density=False)["expr"]
        df["expr"] = (expr_1["expr"] - expr_2["expr"]) / expr_2["expr"]
        df["expr_abs"] = df["expr"].abs()
        df.sort_values("expr_abs", ascending=False, inplace=True)
        df.drop("expr_abs", axis=1, inplace=True)
        return df
