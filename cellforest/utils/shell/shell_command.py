import gc

import logging
import os
import shlex
from pathlib import Path
from subprocess import check_call
from typing import Union

from cellforest import get_current_config


def process_shell_command(command_string: str, logs_dir: Union[str, Path], logfile_prefix: str):
    os.makedirs(logs_dir, exist_ok=True)
    out_path = os.path.join(logs_dir, f"{logfile_prefix}.out")
    err_path = os.path.join(logs_dir, f"{logfile_prefix}.err")
    logger = logging.getLogger(f"{logfile_prefix}:shell command")
    logger.info(f"Starting external process {logfile_prefix}")
    logger.info(f"Running command: {command_string}")
    logger.info(f"STDOUT -> {out_path}")
    logger.info(f"STDERR -> {err_path}")
    config = get_current_config()
    shell_mode = config.get("shell_mode", "subprocess")
    if shell_mode == "subprocess":
        check_call(shlex.split(command_string), stdout=open(out_path, "w"), stderr=open(err_path, "w"))
    elif shell_mode == "os":
        command_string += f" > {out_path} 2> {err_path}"
        ret_code = os.system(command_string)
        if ret_code != 0:
            raise ChildProcessError(f"Error: ret_code {ret_code} raised by command_string: {command_string}")
    logger.info(f"Finished")
    gc.collect()


def shell_command(command_string: str):
    check_call(shlex.split(command_string))
