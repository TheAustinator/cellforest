from functools import wraps
from typing import Union

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix

from cellforest.structures.exceptions import CellsNotFound, GenesNotFound


class Counts(csr_matrix):
    FEATURES_COLUMNS = ["ensgs", "genes"]
    SUPER_METHODS = [
        "arcsin",
        "arcsinh",
        "arctan",
        "arctanh",
        "asformat",
        "asfptype",
        "astype",
        "ceil",
        "conj",
        "conjugate",
        "copy",
        "deg2rad",
        "eliminate_zeros",
        "expm1",
        "floor",
        "log1p",
        "maximum",
        "minimum",
        "multiply",
        "power",
        "rad2deg",
        "rint",
        "sign",
        "sin",
        "sinh",
        "sqrt",
        "tan",
        "tanh",
        "trunc",
    ]

    def __init__(self, matrix, cell_ids, features, **kwargs):
        # TODO: make a get_counts function that just takes the directory
        super().__init__(matrix, **kwargs)
        self.matrix = matrix
        self.features = features.iloc[:, :2].copy()
        self.features.columns = self.FEATURES_COLUMNS
        self._idx = self._convert_to_series(cell_ids)
        self._ids = self._index_col_swap(cell_ids.copy())

    @property
    def genes(self):
        return self.features["genes"]

    @property
    def ensgs(self):
        return self.features["ensgs"]

    @property
    def index(self):
        return self._idx

    @property
    def columns(self):
        return self.genes

    @property
    def cell_ids(self):
        return self.index

    @property
    def _genes_names(self):
        return self._index_col_swap(self.genes.copy(), "genes")

    @property
    def _ensgs_names(self):
        return self._index_col_swap(self.ensgs.copy(), "ensgs")

    def vstack(self, other):
        raise NotImplementedError()

    def hstack(self, other):
        raise NotImplementedError()

    def to_df(self):
        return pd.DataFrame(self.todense(), columns=self.columns, index=self.index)

    def dropna(self, axis=None):
        if axis is None:
            return self.dropna(axis=0).dropna(axis=1)
        sum_axis = int(not bool(axis))
        selector = np.asarray(self.sum(axis=sum_axis)).flatten().astype(bool)
        if axis == 0:
            return self[selector]
        else:
            return self[:, selector]

    @classmethod
    def to_cellranger(cls, output_dir):
        """Save in 10X Cellranger output format"""
        raise NotImplementedError()

    @classmethod
    def from_cellranger(cls, cellranger_dir):
        """Load from 10X Cellranger output format"""
        raise NotImplementedError()

    @classmethod
    def load(cls, filepath):
        """Load from pickle"""
        raise NotImplementedError()

    @classmethod
    def save(cls, filepath):
        """Save as pickle"""
        raise NotImplementedError()

    def __getitem__(self, key):
        if isinstance(key, tuple):
            return self._2d_slice(key)
        else:
            return self._cell_slice(key)

    def _2d_slice(self, key):
        """Slice rows and columns (cells and genes)"""
        gene_sliced = self._gene_slice(key[1])
        cell_sliced = gene_sliced[key[0]]
        return cell_sliced

    def _cell_slice(self, key):
        """Slice rows (cells)"""
        try:
            key = self._convert_key(key, self._ids)
        except KeyError:
            raise CellsNotFound(self._ids, key)
        if isinstance(key, slice):
            cell_ids = pd.DataFrame(self._idx[key]).reset_index(drop=True)
        else:
            cell_ids = pd.DataFrame(self._idx.reindex(key)).reset_index(drop=True)
        mat = csr_matrix(self.matrix)[key]
        return self.__class__(mat, cell_ids, self.features)

    def _gene_slice(self, key):
        """Slice columns (genes) with either gene names or ensemble names"""
        key = self._genes_convert_key(key)
        mat = csr_matrix(self.matrix)[:, key]
        if isinstance(key, slice):
            features = self.features[key]
        else:
            genes = pd.DataFrame(self.genes.reindex(key)).reset_index(drop=True)
            ensgs = pd.DataFrame(self.ensgs.reindex(key)).reset_index(drop=True)
            if len(ensgs) > len(genes):
                features = self.features[self.features.ensgs.isin(key)]
            else:
                features = self.features[self.features.genes.isin(key)]
        return self.__class__(mat, self._idx, features)

    def _genes_convert_key(self, key):
        """"""
        try:
            key = self._convert_key(key, self._genes_names)
        except KeyError:
            try:
                key = self._convert_key(key, self._ensgs_names)
            except KeyError:
                genes_err = GenesNotFound(self._genes_names, key)
                ensgs_err = GenesNotFound(self._ensgs_names, key)
                if len(ensgs_err.missing) < len(genes_err.missing):
                    raise ensgs_err
                else:
                    raise genes_err
        return key

    @staticmethod
    def _convert_key(key, df):
        """Slice index dataframe with key and convert to integer indices"""
        if isinstance(key, (pd.Series, pd.Index, np.ndarray)):
            key = key.tolist()
        if isinstance(key, list):
            if isinstance(key[0], str):
                # gene names are duplicated, ensgs aren't
                if df.index.duplicated().any():
                    df_temp = df.reset_index()
                    key_rows = df_temp[df_temp[df_temp.columns[0]].isin(key)]
                else:
                    key_rows = df.reindex(key)
                key = key_rows.dropna()["i"].astype(int).tolist()
        elif isinstance(key, str):
            key = [df.loc[key]["i"].tolist()]
        elif isinstance(key, int):
            key = [key]
        else:
            return key
        if len(key) == 0:
            raise KeyError("No matching indices")
        return key

    @staticmethod
    def _check_key(key, df):
        """
        NOTE: not currently used because metedata includes cells that were filtered out by cellranger
        """
        intersection = len(set(key).intersection(set(df.index.tolist()))) / len(key)
        if intersection < 1:
            import ipdb

            ipdb.set_trace()
            raise KeyError(f"some of provided keys missing from counts matrix. Intersection: {intersection}")

    @staticmethod
    def _index_col_swap(df, col: Union[str, int] = 0, new_index_colname="i"):
        """Swaps column with index of DataFrame"""
        df = df.copy()
        if isinstance(df, pd.Series):
            df = pd.DataFrame(df)
        df[new_index_colname] = df.index
        df.index = df[col]
        df.drop(columns=col, inplace=True)
        return df

    @staticmethod
    def _convert_to_series(df):
        """If a dataframe, convert to series"""
        if isinstance(df, pd.DataFrame):
            df = df.iloc[:, 0].copy()
        elif not isinstance(df, pd.Series):
            raise TypeError(f"Must be dataframe not series {type(df)}")
        return df

    def __repr__(self):
        return f"{self.__class__}: [cell_ids x genes] matrix\n" + csr_matrix.__repr__(self)

    @staticmethod
    def wrap_super(func):
        """Wrapper to pass scipy matrix methods through to .matrix attribute"""

        @wraps(func)
        def wrapper(counts, *args, **kwargs):
            matrix = func(counts.matrix, *args, **kwargs)
            return counts.__class__(matrix, counts.cell_ids, counts.genes)

        return wrapper

    @staticmethod
    def decorate(method_names):
        """
        Wrap a list of scipy matrix `method_names` with `wrap_super` and
        re-tether them to class
        """
        for name in method_names:
            super_method = getattr(csr_matrix, name)
            wrapped_method = Counts.wrap_super(super_method)
            setattr(Counts, name, wrapped_method)


Counts.decorate(Counts.SUPER_METHODS)
