import shutil

from dataforest.hooks import hook


@hook
def hook_clean_unversioned(dp):
    if dp.branch.unversioned:
        dp.logger.info(f"Removing output files for {dp.process_name} name due to unversioned CellBranch")
        shutil.rmtree(str(dp.branch[dp.process_name].path))
