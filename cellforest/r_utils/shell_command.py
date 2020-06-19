import logging
import os
import shlex
from subprocess import check_call


def shell_command(command_string, working_dir, process_name):
    out_path = os.path.join(working_dir, f"{process_name}.out")
    err_path = os.path.join(working_dir, f"{process_name}.err")

    logger = logging.getLogger(f"{process_name}:shell command")
    logger.info(f"Starting external process {process_name}")
    logger.info(f"Running command: {command_string}")
    logger.info(f"STDOUT -> {out_path}")
    logger.info(f"STDERR -> {err_path}")

    check_call(shlex.split(command_string), stdout=open(out_path, "w"), stderr=open(err_path, "w"))
    logger.info(f"Finished")
