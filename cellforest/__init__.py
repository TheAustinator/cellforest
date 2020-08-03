from cellforest.structures.counts.Counts import Counts
from .templates.CellBranch import CellBranch
from .templates.SpecSC import SpecSC as Spec

from dataforest import get_current_config, update_config

from_metadata = CellBranch.from_metadata
from_input_dirs = CellBranch.from_input_dirs
load = CellBranch.load
