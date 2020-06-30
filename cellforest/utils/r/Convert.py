from pathlib import Path

from cellforest.utils import r
from cellforest.utils.r.shell_command import shell_command

PICKLE_TO_RDS_SCRIPT = Path(r.__file__).parent / "pickle_to_rds.R"


class Convert:
    @staticmethod
    def pickle_to_rds(pickle_path, meta_path, output_rds_path):
        arg_list = [pickle_path, meta_path, output_rds_path]
        Convert._run_r_script(PICKLE_TO_RDS_SCRIPT, arg_list)

    @staticmethod
    def pickle_to_rds_dir(file_dir):
        file_dir = Path(file_dir)
        pickle_path = file_dir / "rna.pickle"
        meta_path = file_dir / "meta.tsv"
        output_rds_path = file_dir / "rna.rds"
        Convert.pickle_to_rds(pickle_path, meta_path, output_rds_path)

    @staticmethod
    def rds_to_pickle():
        raise NotImplementedError()

    @staticmethod
    def _run_r_script(script_path: str, arg_list: list):
        command_string = f"Rscript {str(script_path)} {' '.join(map(str, arg_list))}"
        shell_command(command_string=command_string)
