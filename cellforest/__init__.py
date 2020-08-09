from pathlib import Path

from dataforest import get_current_config, update_config

from cellforest.structures.counts.Counts import Counts
from cellforest.templates.CellBranch import CellBranch
from cellforest.templates.CellInterface import CellInterface
from cellforest.templates.CellTree import CellTree


from_sample_metadata = CellInterface.from_sample_metadata
from_input_dirs = CellInterface.from_input_dirs
load = CellInterface.load

update_config(Path(__file__).parent / "config/default_config.yaml")
