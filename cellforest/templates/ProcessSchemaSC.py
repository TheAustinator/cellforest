import json
from operator import ge, le

from dataforest.templates.ProcessSchema import ProcessSchema
from pathlib import Path

from cellforest.templates.ReaderMethodsSC import ReaderMethodsSC
from cellforest.templates.WriterMethodsSC import WriterMethodsSC


class ProcessSchemaSC(ProcessSchema):
    CONFIG = json.load(Path(__file__).parent / "config/process_schema.json")
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
    FILE_MAP = CONFIG["file_map"]
    READER_METHODS = ReaderMethodsSC
    WRITER_METHODS = WriterMethodsSC
    R_SCRIPTS_DIR = Path(__file__).parent.absolute() / "r_scripts"
    R_FILENAMES = CONFIG["r_filenames"]
    R_FILEPATHS = {k: R_SCRIPTS_DIR / v for k, v in R_FILENAMES.items()}
    TEMP_METADATA_FILENAME = "temp_cell_metadata.tsv"


def get_latest_analysis_run(date):
    pass


def get_analysis_runs_prior_to_lane(lane):
    pass
