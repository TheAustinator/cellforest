import json
from operator import ge, le

from dataforest.templates.ProcessSchema import ProcessSchema
from pathlib import Path

from cellforest.templates.ReaderMethodsSC import ReaderMethodsSC
from cellforest.templates.WriterMethodsSC import WriterMethodsSC


class ProcessSchemaSC(ProcessSchema):
    with open(Path(__file__).parent.parent / "config/process_schema.json") as f:
        CONFIG = json.load(f)
    # TODO: check that process names in PROCESS_HIERARCHY, PARAM_NAMES, and FILE_MAP are consistent
    # fmt: off

    PROCESS_HIERARCHY = CONFIG["process_hierarchy"]
    # fmt: on

    PARAM_NAMES = CONFIG["param_names"]

    PROCESS_NAMES = list(PARAM_NAMES.keys())
    SUBSET_PROXIES = {
        "max_lane_id": (le, "lane_id"),
        "min_lane_id": (ge, "lane_id"),
    }
    STANDARD_FILES = CONFIG["standard_files"]
    FILE_MAP = CONFIG["file_map"]
    READER_METHODS = ReaderMethodsSC
    WRITER_METHODS = WriterMethodsSC
    R_SCRIPTS_DIR = Path(__file__).parent.parent.absolute() / "modules"
    R_FILENAMES = CONFIG["r_filenames"]

    # noinspection
    def _get_r_filepaths(scripts_dir, r_filenames):
        # has to be a function due to scoping issue with class var dict comprehension
        return {k: scripts_dir / v for k, v in r_filenames.items()}

    R_FILEPATHS = _get_r_filepaths(R_SCRIPTS_DIR, R_FILENAMES)
    TEMP_METADATA_FILENAME = "temp_cell_metadata.tsv"


def get_latest_analysis_run(date):
    pass


def get_analysis_runs_prior_to_lane(lane):
    pass
