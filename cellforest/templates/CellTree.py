from dataforest.core.DataTree import DataTree

from cellforest import CellBranch
from cellforest.structures.CellDistributedContainer import CellDistributedContainer
from cellforest.templates.CellBase import CellBase


class CellTree(CellBase, DataTree):
    _BRANCH_CLASS = CellBranch

    @property
    def rna(self):
        # TODO: refactor to soft-coded on `DataTree` using `ASSAY_OPTIONS` attr
        return CellDistributedContainer([branch.rna for branch in self.branches])
