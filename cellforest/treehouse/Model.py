from dataforest.treehouse.model import DataFrameModel


class CellModel(DataFrameModel):
    def __init__(self, tree):
        self.tree = tree
        super().__init__(self.tree.meta)

    @property
    def branch(self):
        return self.tree.branches[0]

    @property
    def df(self):
        return self.branch.meta

    @property
    def rna(self):
        return self.branch.rna
