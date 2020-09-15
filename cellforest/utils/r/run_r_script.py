from cellforest.utils.shell.shell_command import process_shell_command, shell_command


def run_process_r_script(branch: "CellBranch", r_script_filepath: str, arg_list: list, logfile_prefix: str):
    """
    Runs an R script for a process, which additionally entails outputting log
    files
    """
    command_string = f"Rscript {r_script_filepath} {' '.join(map(str, arg_list))}"
    logs_dir = str(branch[logfile_prefix].logs_path)
    process_shell_command(command_string, logs_dir, logfile_prefix)


def run_r_script(script_path: str, arg_list: list):
    command_string = f"Rscript {str(script_path)} {' '.join(map(str, arg_list))}"
    shell_command(command_string)
