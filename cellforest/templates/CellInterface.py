from dataforest.core.Interface import Interface

from cellforest.templates.CellBranch import CellBranch
from cellforest.templates.CellTree import CellTree


class CellInterface(Interface):
    BRANCH_CLASS = CellBranch
    TREE_CLASS = CellTree
