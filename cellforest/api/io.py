from lazy_import import lazy_module, lazy_callable
# from anndata._io.specs.registry import read_elem
read_elem = lazy_callable("anndata._io.specs.registry.read_elem")
from h5py import File
import pandas as pd
from scanpy import read_10x_h5, read_10x_mtx


def read_10x(path):
    try:
        ad = read_10x_h5(path)
    except OSError:
        try:
            ad = read_10x_mtx(str(path)[:-3])
        except OSError:
            raise OSError(f"Neither .h5 or .mtx found for {path}")
    return ad


def read_obs(h5ad_path) -> pd.DataFrame:
    with File(h5ad_path, "r") as f:
        return read_elem(f["obs"])
