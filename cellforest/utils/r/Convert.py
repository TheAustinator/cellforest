import sys

from pathlib import Path

from cellforest.structures.counts import build_counts_store
from cellforest.utils import r
from cellforest.utils.shell.shell_command import shell_command

_PYTHON_EXECUTABLE = sys.executable
_R_UTILS_DIR = Path(r.__file__).parent
_PICKLE_TO_RDS_SCRIPT = _R_UTILS_DIR / "pickle_to_rds.R"
_RDS_TO_PICKLE_SCRIPT = _R_UTILS_DIR / "rds_to_pickle.R"
_COUNTS_STORE_MODULE = build_counts_store.__file__


class Convert:
    # TODO: rewrite using stem rather than hard-coding filenames (used in Counts, too)
    @staticmethod
    def pickle_to_rds(pickle_path, meta_path, output_rds_path):
        arg_list = [pickle_path, meta_path, output_rds_path, _PYTHON_EXECUTABLE]
        Convert._run_r_script(_PICKLE_TO_RDS_SCRIPT, arg_list)

    @staticmethod
    def pickle_to_rds_dir(file_dir):
        pickle_path, meta_path, output_rds_path = Convert._get_std_paths(file_dir)
        Convert.pickle_to_rds(pickle_path, meta_path, output_rds_path)

    @staticmethod
    def rds_to_pickle(rds_path, output_meta_path, output_pickle_path):
        arg_list = [rds_path, output_meta_path, output_pickle_path, _COUNTS_STORE_MODULE, _PYTHON_EXECUTABLE]
        Convert._run_r_script(_RDS_TO_PICKLE_SCRIPT, arg_list)

    @staticmethod
    def rds_to_pickle_dir(file_dir):
        output_pickle_path, output_meta_path, rds_path = Convert._get_std_paths(file_dir)
        Convert.rds_to_pickle(rds_path, output_meta_path, output_pickle_path)

    @staticmethod
    def _run_r_script(script_path: str, arg_list: list):
        command_string = f"Rscript {str(script_path)} {' '.join(map(str, arg_list))}"
        shell_command(command_string=command_string)

    @staticmethod
    def _get_std_paths(file_dir):
        file_dir = Path(file_dir)
        pickle_path = file_dir / "rna.pickle"
        meta_path = file_dir / "meta.tsv"
        rds_path = file_dir / "rna.rds"
        return pickle_path, meta_path, rds_path
