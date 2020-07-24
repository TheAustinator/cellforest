import os
import pickle
from pathlib import Path

from cellforest.structures.counts.CountsStore import CountsStore


def build_counts_store(matrix, cell_ids, features, save_path=None):
    """
    Method for saving counts matrix data as a pickle file. Intended for use in
    python or R (via Reticulate).
    Args:
        matrix: sparse matrix
        cell_ids:
        features:
        save_path:

    Returns:

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
