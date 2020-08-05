from dataforest.core.DataTree import DataTree

from cellforest import CellBranch
from cellforest.templates.CellDataBase import CellDataBase


class CellTree(CellDataBase, DataTree):
    BRANCH_CLASS = CellBranch
