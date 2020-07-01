from pathlib import Path

from cellforest.utils import r
from cellforest.utils.r.shell_command import shell_command

r_utils = Path(r.__file__).parent
_PICKLE_TO_RDS_SCRIPT = r_utils / "pickle_to_rds.R"
_RDS_TO_PICKLE_SCRIPT = r_utils / "rds_to_pickle.R"
_COUNTS_STORE_MODULE = r_utils.parent.parent / "structures/build_counts_store.py"


class Convert:
    @staticmethod
    def pickle_to_rds(pickle_path, meta_path, output_rds_path):
        arg_list = [pickle_path, meta_path, output_rds_path]
        Convert._run_r_script(_PICKLE_TO_RDS_SCRIPT, arg_list)

    @staticmethod
    def pickle_to_rds_dir(file_dir):
        pickle_path, meta_path, output_rds_path = Convert._get_std_paths(file_dir)
        Convert.pickle_to_rds(pickle_path, meta_path, output_rds_path)

    @staticmethod
    def rds_to_pickle(rds_path, output_meta_path, output_pickle_path):
        arg_list = [rds_path, output_meta_path, output_pickle_path, _COUNTS_STORE_MODULE]
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
