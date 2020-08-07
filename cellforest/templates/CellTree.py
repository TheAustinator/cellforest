from dataforest.core.DataTree import DataTree

from cellforest import CellBranch
from cellforest.templates.CellBase import CellBase


class CellTree(CellBase, DataTree):
    BRANCH_CLASS = CellBranch
