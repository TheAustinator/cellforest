import shutil

from dataforest.hooks import hook


@hook
def hook_clean_unversioned(dp):
    if dp.forest.unversioned:
        dp.logger.info(f"Removing output files for {dp.process_name} name due to unversioned CellForest")
        shutil.rmtree(str(dp.forest[dp.process_name].path))
