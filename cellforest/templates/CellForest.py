from copy import deepcopy
from pathlib import Path
from typing import Optional, Union, List

from dataforest.core.DataForest import DataForest
from dataforest.utils.utils import label_df_partitions, update_recursive
import pandas as pd

from cellforest.structures.Counts import Counts
from cellforest.templates.PlotMethodsSC import PlotMethodsSC
from cellforest.templates.ProcessMethodsSC import ProcessMethodsSC
from cellforest.templates.ProcessSchemaSC import ProcessSchemaSC
from cellforest.templates.ReaderMethodsSC import ReaderMethodsSC
from cellforest.templates.SpecSC import SpecSC
from cellforest.templates.WriterMethodsSC import WriterMethodsSC


class CellForest(DataForest):
    """
    DataForest for scRNAseq processed data. The `process_hierarchy` currently starts
    at `combine`, where non-normalized counts data is combined.

    A path through specific `process_runs` of processes in the
    `process_hierarchy` are specified in the `spec_dict`, according to the
    specifications of `dataforest.Spec`. Any root level (not under a process
    name in `spec_dict`) `subset`s or `filter`s are applied to `counts` and
    `meta`, which are the preferred methods for accessing cell metadata and
    the normalized counts matrix
    """

    ROOT_LEVEL_COMPARTMENTS = {
        "subset",
    }
    PLOT_METHODS = PlotMethodsSC
    PROCESS_METHODS = ProcessMethodsSC
    SPEC_CLASS = SpecSC
    SCHEMA_CLASS = ProcessSchemaSC
    READER_METHODS = ReaderMethodsSC
    WRITER_METHODS = WriterMethodsSC
    READER_KWARGS_MAP = {
        "dim_reduce": {
            "pca_embeddings": {"header": "infer"},
            "pca_loadings": {"header": "infer"},
            "umap_embeddings": {"header": "infer", "index_col": 0},
        },
        "combine": {"cell_metadata": {"header": 0}},
        # 'normalize': {
        #     'cell_ids': {'index_col': 0}},
        "cluster": {"clusters": {"index_col": 0}},
        "diffexp": {"diffexp_result": {"header": 0}},
    }
    METADATA_NAME = "meta"
    COPY_KWARGS = {**DataForest.COPY_KWARGS, "unversioned": "unversioned"}

    def __init__(self, root_dir, spec_dict=None, verbose=False, meta=None, unversioned=None):
        super().__init__(root_dir, spec_dict, verbose)
        self._rna = None
        self._meta_unfiltered = None
        if meta is not None:
            meta = meta.copy()
        self._meta = self.get_cell_meta(meta)
        # TODO: use this to augment strings of output directories so manual tinkers don't
        #   affect downstream processing
        if meta is not None and unversioned is None:
            self._unversioned = True
        else:
            self._unversioned = bool(unversioned)
        if self.unversioned:
            self.logger.warning(f"Unversioned DataForest")

    @property
    def samples(self) -> pd.DataFrame:
        """
        Hierarchical categorization of all samples in dataset with cell counts.
        The canonical use case would be to use it on a broad DataForest to choose
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
        if self._rna is None:
            # TODO: set to use normalized if exists by default -- see old version
            # TODO: soft code filenames
            counts_path = self.root_dir / "counts.pickle"
            if not counts_path.exists:
                raise FileNotFoundError(
                    f"Ensure that you initialized the root directory with CellForest.from_metadata or "
                    f"CellForest.from_input_dirs. Not found: {counts_path}"
                )
            self._rna = Counts.load(counts_path)
            self._rna = self._rna[self.meta.index]
        return self._rna

    @property
    def meta(self) -> pd.DataFrame:
        """
        Interface for cell metadata, which is derived from the sample
        metadata and the scrnaseq experimental data. Available UMAP embeddings
        and cluster identifiers will be included, and the data will be subset,
        filtered, and partitioned based on the specifications in `self.spec`.
        Primarily for this reason, this is the preferred interface to metadata
        over direct file access.
        """
        # TODO: add embeddings and cluster ids
        if self._meta is None:
            self._meta = self.get_cell_meta()
        elif "dim_reduce" in self.spec and "UMAP_1" not in self._meta.columns:
            if self["dim_reduce"].done:
                self._meta = self.get_cell_meta()
        elif "cluster" in self.spec and "cluster_id" not in self._meta.columns:
            if self["cluster"].done:
                self._meta = self.get_cell_meta()
        return self._meta

    @property
    def unversioned(self):
        return self._unversioned

    @property
    def meta_unfiltered(self) -> pd.DataFrame:
        # TODO: not used anywhere, figure out use and add docstring or delete
        return self._meta_unfiltered

    def set_partition(self, process_name: Optional[str] = None, encodings=True):
        """Add columns to metadata to indicate partition from spec"""
        columns = self.spec[process_name]["partition"]
        self._meta = label_df_partitions(self.meta, columns, encodings)

    def groupby(self, by: Union[str, list, set, tuple], **kwargs):
        """
        Operates like a pandas groupby, but does not return a GroupBy object,
        and yields (name, DataForest), where each DataForest is subset according to `by`,
        which corresponds to columns of `self.meta`.
        This is useful for batching analysis across various conditions, where
        each run requires an DataForest.
        Args:
            by: variables over which to group (like pandas)
            **kwargs: for pandas groupby on `self.meta`

        Yields:
            name: values for DataForest `subset` according to keys specified in `by`
            forest: new DataForest which inherits `self.spec` with additional `subset`s
                from `by`
        """
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
            # forest._meta = df
            yield name, forest

    def copy(self, reset: bool = False, **kwargs):
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

    def _get_compartment_updated(self, compartment_name: str, update: dict) -> "CellForest":
        """

        """
        if compartment_name in self.ROOT_LEVEL_COMPARTMENTS:
            spec = update_recursive(self.spec, update, inplace=False)
        else:
            spec = update_recursive(self.spec, {compartment_name: update}, inplace=False)
        forest = self.copy(spec_dict=spec)
        bool_selector = pd.concat([forest.meta[key] == value for key, value in update.items()], axis=1).all(axis=1)
        if sum(bool_selector) == 0:
            import ipdb

            ipdb.set_trace()
        if compartment_name == "subset":
            forest._meta = forest.meta[bool_selector]
        elif compartment_name == "filter":
            forest._meta = forest.meta[~bool_selector]
        else:
            raise ValueError()
        forest._rna = forest.counts[forest.meta.index]
        return forest

    def get_cell_meta(self, df=None):
        # TODO: MEMORY DUPLICATION - we want to keep file access pure?
        if df is None:
            # TODO: fix this
            try:
                # df = self.f["cell_metadata"].copy()
                df = pd.read_csv(self.root_dir / "meta.tsv", sep="\t", index_col=0)
            except FileNotFoundError:
                df = pd.DataFrame(self.rna.cell_ids.copy())
                df.columns = ["cell_id"]
                df.index = df["cell_id"]
                df.drop(columns=["cell_id"], inplace=True)
            df.replace(" ", "_", regex=True, inplace=True)
            if "to_bucket_var" in df and "bucketed_var" not in df:
                df["bucketed_var"] = pd.cut(df["to_bucket_var"], bins=(0, 20, 40, 60, 80), labels=(10, 30, 50, 70),)
            if "str_var_preprocessed" in df and "str_var_processed" not in df:
                df["str_var_processed"] = df["str_var_preprocessed"].str.extract(r"([A-Z]\d)")
            # TODO: fill in once `process_run.done` feature is ready
            df = self._meta_add_downstream_data(df)
        df = self._subset_filter(df, self.spec, self.schema)
        if self.spec.partition_set:
            df = label_df_partitions(df, self.spec.partition_set, encodings=True)
        return df

    def _meta_add_downstream_data(self, df):
        done = set()
        if "cluster" in self.spec:
            if self["cluster"].done:
                done.update({"normalize", "dim_reduce", "cluster"})
        if not done and "dim_reduce" in self.spec:
            if self["dim_reduce"].done:
                done.update({"normalize", "dim_reduce"})
        if not done and "normalize" in self.spec:
            if self["normalize"].done:
                done.update({"normalize"})
        if "cluster" in done:
            clusters = self.f["cluster"]["clusters"].copy()
            clusters.rename(columns={1: "cluster_id"}, inplace=True)
            df = df.merge(clusters, how="left", left_index=True, right_index=True)
            df["cluster_id"] = df["cluster_id"].astype(pd.Int16Dtype())
        if "dim_reduce" in done:
            df = df.merge(self.f["dim_reduce"]["umap_embeddings"], how="left", left_index=True, right_index=True)
        if "normalize" in done:
            df = df[df.index.isin(self.f["normalize"]["cell_ids"][0])]
        # TODO: temp during param mismatch
        try:
            df = df[df.index.isin(self.f["normalize"]["cell_ids"][0])]
        except Exception:
            self.logger.info("Could not find filtered cell ids. Using all cells in metadata")
        return df

    @staticmethod
    def _combine_datasets(
        root_dir: Union[str, Path],
        metadata: Optional[Union[str, Path, pd.DataFrame]] = None,
        input_paths: Optional[List[Union[str, Path]]] = None,
        metadata_read_kwargs: Optional[dict] = None,
    ):
        """
        Combine files from multiple cellranger output directories into a single
        `Counts` and save it to `root_dir`. If sample metadata is provided,
        replicate each row corresponding to the number of cells in the sample
        such that the number of rows changes from n_samples to n_cells.
        """
        root_dir = Path(root_dir)
        if (input_paths and metadata) or (input_paths is None and metadata is None):
            raise ValueError("Must specify exactly one of `input_dirs` or `metadata`")
        elif metadata is not None:
            if isinstance(metadata, (str, Path)):
                metadata_read_kwargs = {"sep": "\t"} if not metadata_read_kwargs else metadata_read_kwargs
                metadata = pd.read_csv(metadata, **metadata_read_kwargs)
            if "path" not in metadata.columns:
                raise ValueError("metadata must contain column: `path` with reference to cellranger outputs")
            paths = metadata["path"]
            counts_list = [Counts.from_cellranger(dir_) for dir_ in paths]
            cells_per_matrix = [counts.shape[0] for counts in counts_list]
            meta = metadata.loc[metadata.index.repeat(cells_per_matrix)].reset_index(drop=True)
        else:
            counts_list = [Counts.from_cellranger(dir_) for dir_ in input_paths]
            meta = None
        counts = Counts.concatenate(counts_list)
        counts.save(root_dir / "counts.pickle")
        if meta is not None:
            meta.index = counts.cell_ids
            meta.to_csv(root_dir / "meta.tsv", sep="\t")
        return dict()
