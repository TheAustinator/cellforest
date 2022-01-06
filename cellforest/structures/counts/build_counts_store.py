import os
import pickle
from pathlib import Path
from typing import Union, Optional

from cellforest.structures.counts.CountsStore import CountsStore
import numpy as np
import pandas as pd
from scipy.sparse.base import spmatrix


def build_counts_store(
    matrix: Union[np.ndarray, spmatrix],
    cell_ids: pd.Series,
    features: pd.DataFrame,
    save_path: Optional[Union[str, Path]] = None,
) -> CountsStore:
    """
    Method for saving counts matrix data as a pickle file. Intended for use in
    python or R (via Reticulate).
    Args:
        matrix: see `cellforest.structures.counts.Counts.Counts`
        cell_ids: see `cellforest.structures.counts.Counts.Counts`
        features: see `cellforest.structures.counts.Counts.Counts`
        save_path: path to store pickle object

    Returns:
        store: dummy object to hold data for pickling
    """
    store = CountsStore()
    store.matrix = matrix
    store.cell_ids = cell_ids
    store.features = features
    if save_path:
        save_path = Path(save_path)
        os.makedirs(str(save_path.parent), exist_ok=True)
        with open(str(save_path), "wb") as f:
            pickle.dump(store, f, protocol=pickle.HIGHEST_PROTOCOL)
    return store
