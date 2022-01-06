import pandas as pd


class SliceErrorBase(KeyError):
    INDEX_TYPE = "index items"

    def __init__(self, index, key):
        key = pd.Series(key)
        self.missing = key[~key.isin(list(index))].tolist()

    def __str__(self):
        return f"The following {self.INDEX_TYPE}s were not found: {self.missing}"


class CellsNotFound(SliceErrorBase):
    INDEX_TYPE = "`cell_ids`"


class GenesNotFound(SliceErrorBase):
    INDEX_TYPE = "`genes`"
