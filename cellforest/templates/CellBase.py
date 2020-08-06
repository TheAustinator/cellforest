import os
from pathlib import Path
from typing import Union, Optional, List

import pandas as pd
from dataforest.core.DataBase import DataBase

from cellforest.utils.cellranger.DataMerge import DataMerge


class CellBase(DataBase):
    _ASSAY_OPTIONS = ["rna", "vdj", "surface", "antigen", "cnv", "atac", "spatial", "crispr"]
    _DEFAULT_CONFIG = Path(__file__).parent.parent / "config/default_config.yaml"

    @staticmethod
    def _combine_datasets(
        root: Union[str, Path],
        metadata: Optional[Union[str, Path, pd.DataFrame]] = None,
        input_paths: Optional[List[Union[str, Path]]] = None,
        metadata_read_kwargs: Optional[dict] = None,
        mode: Optional[str] = None,
    ):
        """
        Combine files from multiple cellranger output directories into a single
        `Counts` and save it to `root`. If sample metadata is provided,
        replicate each row corresponding to the number of cells in the sample
        such that the number of rows changes from n_samples to n_cells.
        """
        root = Path(root)
        mode = mode if mode else "rna"
        if (input_paths and metadata) or (input_paths is None and metadata is None):
            raise ValueError("Must specify exactly one of `input_dirs` or `metadata`")
        elif metadata is not None:
            if isinstance(metadata, (str, Path)):
                metadata_read_kwargs = {"sep": "\t"} if not metadata_read_kwargs else metadata_read_kwargs
                metadata = pd.read_csv(metadata, **metadata_read_kwargs)
            prefix = "path_"
            assays = [x[len(prefix) :] for x in metadata.columns if x.startswith(prefix)]
            if len(assays) == 0:
                raise ValueError(
                    f"metadata must contain at least once column named with the prefix, `path_`, and one of the "
                    f"following assays as a suffix: {CellBase._ASSAY_OPTIONS}"
                )
            for assay in assays:
                paths = metadata[f"{prefix}{assay}"].tolist()
                DataMerge.merge_assay(paths, assay, metadata, save_dir=root)
        else:
            DataMerge.merge_assay(input_paths, mode, save_dir=root)
        return dict()

    @staticmethod
    def _get_assays(path):
        # TODO: will have to change once decoupled from pickle (e.g. rds, anndata)
        files = list(filter(lambda x: x.endswith(".pickle"), os.listdir(path)))
        return set(map(lambda x: x.split(".")[0], files))
