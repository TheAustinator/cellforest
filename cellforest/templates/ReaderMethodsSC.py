import pickle

from dataforest.utils.decorators import default_kwargs
import pandas as pd
from pathlib import Path
from scipy import io


class ReaderMethodsSC:
    TSV_DEFAULTS = {"sep": "\t", "header": None}
    TSV_GZ_DEFAULTS = {**TSV_DEFAULTS, "filesystem": "gzip"}
    PICKLE_DEFAULTS = {}

    @staticmethod
    @default_kwargs(TSV_DEFAULTS)
    def tsv(filepath, **kwargs):
        df = pd.read_csv(filepath, **kwargs)
        return df

    @staticmethod
    @default_kwargs(TSV_GZ_DEFAULTS)
    def tsv_gz(filepath, **kwargs):
        # TODO: should be able to pass kwargs directly now, given decorator
        defaults = {"sep": "\t", "filesystem": "gzip", "header": None}
        defaults.update(kwargs)
        kwargs = defaults
        df = pd.read_csv(filepath, **kwargs)
        return df

    @staticmethod
    @default_kwargs(PICKLE_DEFAULTS)
    def pickle(filepath, **kwargs):
        with Path(filepath).open("rb") as f:
            mat = pickle.load(f, **kwargs)
        return mat

    @staticmethod
    def rds(filepath):
        raise NotImplementedError()

    @staticmethod
    def mtx_gz(filepath):
        return io.mmread(str(filepath)).T
