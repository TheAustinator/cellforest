import os
import logging
from pathlib import Path

import pandas as pd
from scipy import io
from scipy.io import mmwrite

from cellforest.utils import compress


class CellRangerIO:
    _MATRIX_BASENAME = "matrix.mtx"
    _FEATURES_BASENAME = "features.tsv"
    _BARCODES_BASENAME = "barcodes.tsv"
    _ALT_FEATURES_BASENAME = "genes.tsv"
    _READ_TSV_KWARGS = {"sep": "\t", "header": None}
    _WRITE_TSV_KWARGS = {"sep": "\t", "header": None, "index": False}

    def __init__(self, cellranger_dir):
        """Initialize for reading only, not for writing"""
        self.logger = logging.getLogger(self.__class__.__name__)
        self.cellranger_dir = Path(cellranger_dir)
        self.files = os.listdir(cellranger_dir)
        self.matrix_filename = self._get_filename(self.files, self._MATRIX_BASENAME)
        self.barcodes_filename = self._get_filename(self.files, self._BARCODES_BASENAME)
        try:
            self.features_filename = self._get_filename(self.files, self._FEATURES_BASENAME)
        except IndexError:
            self.features_filename = self._get_filename(self.files, self._ALT_FEATURES_BASENAME)
            self.logger.warning(f"Old 10X chemistry detected -- features named: {self.features_filename}")
        self._is_matrix_gz = self._is_gz(self.matrix_filename)
        self._is_features_gz = self._is_gz(self.features_filename)
        self._is_barcodes_gz = self._is_gz(self.barcodes_filename)
        self.read_matrix = self._tether(self.read_matrix, self.cellranger_dir / self.matrix_filename)
        self.read_features = self._tether(self.read_features, self.cellranger_dir / self.features_filename)
        self.read_barcodes = self._tether(self.read_barcodes, self.cellranger_dir / self.barcodes_filename)

    @staticmethod
    def read_barcodes(filepath, **kwargs):
        # TODO: convert these to use default_args decorator
        kwargs = {**CellRangerIO._READ_TSV_KWARGS.copy(), **kwargs}
        if str(filepath).endswith(".gz"):
            kwargs.update({"compression": "gzip"})
        df = pd.read_csv(filepath, **kwargs)
        return df

    @staticmethod
    def read_features(filepath, **kwargs):
        kwargs = {**CellRangerIO._READ_TSV_KWARGS.copy(), **kwargs}
        if str(filepath).endswith(".gz"):
            kwargs.update({"compression": "gzip"})
        df = pd.read_csv(filepath, **kwargs)
        return df

    @staticmethod
    def read_matrix(filepath):
        return io.mmread(str(filepath)).T

    @staticmethod
    def write_barcodes(filepath, df, gz=True):
        df.to_csv(filepath, **CellRangerIO._WRITE_TSV_KWARGS)
        CellRangerIO._gzip_if_needed(filepath, gz)

    @staticmethod
    def write_features(filepath, df, gz=True):
        df.to_csv(filepath, **CellRangerIO._WRITE_TSV_KWARGS)
        CellRangerIO._gzip_if_needed(filepath, gz)

    @staticmethod
    def write_matrix(filepath, matrix, gz=True):
        mmwrite(str(filepath), matrix.T)
        CellRangerIO._gzip_if_needed(filepath, gz)

    @staticmethod
    def _get_filename(files, basename):
        return list(filter(lambda x: x.startswith(basename), files))[0]

    @staticmethod
    def _is_gz(filepath):
        return filepath.endswith(".gz")

    @staticmethod
    def _gzip_if_needed(path, gzip):
        if gzip:
            compress(path, keep_orig=False)

    @staticmethod
    def _tether(method, tether_arg):
        def tethered_method(*args, **kwargs):
            return method(tether_arg, *args, **kwargs)

        return tethered_method
