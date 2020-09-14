from cellforest.utils.shell.shell_command import process_shell_command, shell_command


def run_process_r_script(branch: "CellBranch", r_script_filepath: str, arg_list: list, process_name: str):
    """
    Runs an R script for a process, which additionally entails outputting log
    files
    """
    command_string = f"Rscript {r_script_filepath} {' '.join(map(str, arg_list))}"
    working_dir = str(branch[process_name].logs_path)
    process_shell_command(
        command_string=command_string, working_dir=working_dir, process_name=process_name,
    )


def run_r_script(script_path: str, arg_list: list):
    command_string = f"Rscript {str(script_path)} {' '.join(map(str, arg_list))}"
    shell_command(command_string=command_string)
