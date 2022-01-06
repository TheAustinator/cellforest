import warnings

import pandas as pd
from numba import NumbaPerformanceWarning
import umap


def run_umap(
    input_file_path: str,
    n_neighbors: int = 10,
    min_dist: float = 0.5,
    n_components: int = 2,
    metric: str = "euclidean",
    seed: int = 42,
):
    warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)
    pca_df = pd.read_csv(input_file_path, sep="\t")
    umap_handle = umap.UMAP(
        n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed, n_components=n_components, metric=metric,
    )
    umap_matrix = umap_handle.fit(pca_df).embedding_
    umap_df = pd.DataFrame(umap_matrix, columns=[f"UMAP_{idx + 1}" for idx in range(umap_matrix.shape[1])],)
    return umap_df
