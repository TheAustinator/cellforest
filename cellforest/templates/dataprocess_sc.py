import logging
import os
import shutil

from dataforest.core.dataprocess import dataprocess

from cellforest.templates.WriterMethodsSC import WriterMethodsSC


# noinspection PyPep8Naming
class dataprocess_sc(dataprocess):
    logger = logging.getLogger("dataprocess_sc")

    def __init__(
        self, requires: str, comparative: bool = False, overwrite: bool = True, temp_meta: bool = True,
    ):
        super().__init__(requires, comparative, overwrite)
        if temp_meta:
            self.setup_hooks.append(self._hook_store_temp_metadata)
            self.clean_hooks.append(self._hook_clean_temp_metadata)
            self.clean_hooks.append(self._hook_clean_unversioned)

    def _hook_store_temp_metadata(self):
        """
        Stores a temporary metadata file in the current process
        """
        self._metadata_filepath = self.forest[self.process_name].path / self.forest.schema.TEMP_METADATA_FILENAME
        WriterMethodsSC.tsv(
            self._metadata_filepath, self.forest[self.process_name].forest.meta, header=True,
        )

    def _hook_clean_temp_metadata(self):
        os.remove(self._metadata_filepath)

    def _hook_clean_unversioned(self):
        if self.forest.unversioned:
            self.logger.info(f"Removing output files for {self.process_name} name due to unversioned CellForest")
            shutil.rmtree(str(self.forest[self.process_name].path))
