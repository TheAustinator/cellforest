import pandas as pd

from dataforest.decorators import default_kwargs


class WriterMethodsSC:
    TSV_DEFAULTS = {"sep": "\t", "header": None}
    TSV_GZ_DEFAULTS = {**TSV_DEFAULTS, "compression": "gzip"}

    @staticmethod
    @default_kwargs(TSV_DEFAULTS)
    def tsv(filepath: str, obj: pd.DataFrame, **kwargs):
        obj.to_csv(filepath, **kwargs)

    @staticmethod
    @default_kwargs(TSV_GZ_DEFAULTS)
    def tsv_gz(filepath, obj, **kwargs):
        WriterMethodsSC.tsv(filepath, obj, **kwargs)

    @staticmethod
    def pickle(filepath, obj, **kwargs):
        raise NotImplementedError()

    @staticmethod
    def rds(filepath, obj, **kwargs):
        raise NotImplementedError()

    @staticmethod
    def mtx_gz(filepath, obj, **kwargs):
        raise NotImplementedError()
