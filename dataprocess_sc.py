import logging
import os
import shutil

from dataforest.dataprocess import dataprocess
from cellforest.WriterMethodsSC import WriterMethodsSC


# noinspection PyPep8Naming
class dataprocess_sc(dataprocess):
    logger = logging.getLogger("dataprocess_sc")

    def __init__(self, requires: str, comparative: bool = False, overwrite: bool = True, temp_meta: bool = True):
        super().__init__(requires, comparative, overwrite)
        if temp_meta:
            self.setup_hooks.append(self._hook_store_temp_metadata)
            self.clean_hooks.append(self._hook_clean_temp_metadata)
            self.clean_hooks.append(self._hook_clean_unversioned)

    def _hook_store_temp_metadata(self):
        self._metadata_filepath = self.orm[self.process_name].path / self.orm.schema.TEMP_METADATA_FILENAME
        WriterMethodsSC.tsv(self._metadata_filepath, self.orm[self.process_name].orm.meta, header=True)

    def _hook_clean_temp_metadata(self):
        os.remove(self._metadata_filepath)

    def _hook_clean_unversioned(self):
        if self.orm.unversioned:
            self.logger.info(f"Removing output files for {self.process_name} name due to unversioned ORM")
            shutil.rmtree(str(self.orm[self.process_name].path))
