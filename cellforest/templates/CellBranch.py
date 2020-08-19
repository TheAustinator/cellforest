import os
from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List, Tuple

from dataforest.core.DataBranch import DataBranch
from dataforest.core.BranchSpec import BranchSpec
from dataforest.utils.utils import label_df_partitions
import pandas as pd

from cellforest.structures.counts.Counts import Counts
from cellforest.templates.CellBase import CellBase
from cellforest.templates.ReaderMethodsSC import ReaderMethodsSC
from cellforest.templates.WriterMethodsSC import WriterMethodsSC
from cellforest.utils.cellranger.DataMerge import DataMerge


class CellBranch(CellBase, DataBranch):
    """
    DataBranch for scRNAseq processed data. The `process_hierarchy` currently
    starts at `combine`, where non-normalized counts data is combined.

    A path through specific `process_runs` of processes in the
    `process_hierarchy` are specified in the `branch_spec`, according to the
    specifications of `dataforest.Spec`. Any root level (not under a process
    name in `branch_spec`) `subset`s or `filter`s are applied to `counts` and
    `meta`, which are the preferred methods for accessing cell metadata and
    the normalized counts matrix
    """

    READER_METHODS = ReaderMethodsSC
    WRITER_METHODS = WriterMethodsSC
    DATA_FILE_ALIASES = {"rna", "vdj", "surface", "antigen", "cnv", "atac", "spatial", "crispr"}
    READER_KWARGS_MAP = {
        "reduce": {
            "pca_embeddings": {"header": "infer"},
            "pca_loadings": {"header": "infer"},
            "umap_embeddings": {"header": "infer", "index_col": 0},
        },
        "combine": {"cell_metadata": {"header": 0}},
        "cluster": {"clusters": {"index_col": 0}},
        "diffexp": {"diffexp": {"header": 0}},
    }
    _METADATA_NAME = "meta"
    _COPY_KWARGS = {**DataBranch._COPY_KWARGS, "unversioned": "unversioned"}

    def __init__(
        self,
        root: Union[str, Path],
        branch_spec: Optional[Union[list, BranchSpec]] = None,
        verbose: bool = False,
        current_process: Optional[str] = None,
        remote_root: Optional[Union[str, Path]] = None,
        unversioned: Optional[bool] = None,
    ):
        super().__init__(root, branch_spec, verbose, current_process, remote_root)
        self.assays = set()
        self._rna = None
        self._meta = None
        self._unversioned = unversioned

    @property
    def samples(self) -> pd.DataFrame:
        """
        Hierarchical categorization of all samples in dataset with cell counts.
        The canonical use case would be to use it on a broad DataBranch to choose
        a dataset.
        Returns:

        """
        raise NotImplementedError()

    @property
    def rna(self) -> Counts:
        """
        Interface for normalized counts data. It uses the `Counts` wrapper
        around `scipy.sparse.csr_matrix`, which allows for slicing with
        `cell_id`s and `gene_name`s.
        """
        if self._rna is None or not self._rna.index.equals(self.meta.index):
            if self.current_process is not None:
                path_map = self[self.current_process].path_map
                counts_path = path_map["rna"]
            else:
                # TODO: fix hardcoding
                counts_path = self.root / "rna.pickle"
            if not counts_path.exists():
                raise FileNotFoundError(
                    f"Ensure that you initialized the root directory with CellBranch.from_sample_metadata or "
                    f"CellBranch.from_input_dirs. Not found: {counts_path}"
                )
            self._rna = Counts.load(counts_path)
        if not self._rna.index.equals(self.meta.index):
            self._rna = self._rna[self.meta.index]
        return self._rna

    @property
    def vdj(self):
        raise NotImplementedError()

    @property
    def surface(self):
        raise NotImplementedError()

    @property
    def antigen(self):
        raise NotImplementedError()

    @property
    def cnv(self):
        raise NotImplementedError()

    @property
    def atac(self):
        raise NotImplementedError()

    @property
    def spatial(self):
        raise NotImplementedError()

    @property
    def crispr(self):
        raise NotImplementedError()

    def groupby(self, by: Union[str, list, set, tuple], **kwargs) -> Tuple[str, "CellBranch"]:
        """
        Operates like a pandas groupby, but does not return a GroupBy object,
        and yields (name, DataBranch), where each DataBranch is subset according to `by`,
        which corresponds to columns of `self.meta`.
        This is useful for batching analysis across various conditions, where
        each run requires an DataBranch.
        Args:
            by: variables over which to group (like pandas)
            **kwargs: for pandas groupby on `self.meta`

        Yields:
            name: values for DataBranch `subset` according to keys specified in `by`
            branch: new DataBranch which inherits `self.spec` with additional `subset`s
                from `by`
        """
        raise NotImplementedError("currently not functioning")
        if isinstance(by, (tuple, set)):
            by = list(by)
        for (name, df) in self.meta.groupby(by, **kwargs):
            if isinstance(by, list):
                if isinstance(name, (list, tuple)):
                    subset_dict = dict(zip(by, name))
                else:
                    subset_dict = {by[0]: name}
            else:
                subset_dict = {by: name}
            forest = self.get_subset(subset_dict)
            # branch._meta = df
            yield name, forest

    def copy(self, reset: bool = False, **kwargs) -> "CellBranch":
        if kwargs.get("meta", None) is not None:
            kwargs["unversioned"] = True
        if not kwargs:
            kwargs["meta"] = self._meta  # save compute if no modifications
        base_kwargs = self._get_copy_base_kwargs()
        kwargs = {**base_kwargs, **kwargs}
        kwargs = {k: deepcopy(v) for k, v in kwargs.items()}
        if reset:
            kwargs = base_kwargs
        return self.__class__(**kwargs)

    def set_partition(self, process_name: Optional[str] = None, encodings=True):
        """Add columns to metadata to indicate partition from branch_spec"""
        columns = self.spec[process_name]["_PARTITION_"]
        self._meta = label_df_partitions(self.meta, columns, encodings)

    def _get_meta(self, process_name: str) -> pd.DataFrame:
        """
        Read in cell metadata and performs modifications:
            - replace any spaces with underscores
            - merge metadata from process precursors and current process
            - performs data operations (subset, filter, partition)
        Args:
            df: if provided, skip first two steps and go straight to data ops
        Returns:
            df: modified metadata dataframe
        """
        try:
            df = pd.read_csv(self["root"].path_map["meta"], sep="\t", index_col=0)
        except FileNotFoundError:
            df = pd.DataFrame(self.rna.cell_ids.copy())
            df.columns = ["cell_id"]
            df.index = df["cell_id"]
            df.drop(columns=["cell_id"], inplace=True)
        if process_name is not None:
            precursor_names = self.spec.get_precursors_lookup(incl_current=True)[process_name]
            for precursor_name in precursor_names:
                try:
                    process_meta = self[precursor_name].process_meta
                except FileNotFoundError:
                    pass
                else:
                    intersect_cols = set(df.columns).intersection(set(process_meta.columns))
                    process_meta.drop(intersect_cols, axis=1, inplace=True)
                    df = df.merge(process_meta, left_index=True, right_index=True)
        return df
