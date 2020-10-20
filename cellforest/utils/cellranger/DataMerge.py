import gc
import os

from joblib import Parallel, delayed

from cellforest import Counts


class DataMerge:
    @staticmethod
    def merge_assay(paths, mode, metadata=None, save_dir=None):
        method = getattr(DataMerge, f"_merge_{mode}")
        ret = method(paths, metadata, save_dir)
        gc.collect()
        return ret

    @staticmethod
    def _merge_rna(paths, metadata, save_dir, id_col="lane_id"):
        """"""
        pool = Parallel(n_jobs=-2)
        # TODO: significant memory leakage -- maybe make an optional kwarg
        rna_list = pool(delayed(Counts.from_cellranger)(path) for path in paths)
        widths = list(map(lambda x: x.shape[1], rna_list))
        if len(set(widths)) > 1:
            raise ValueError(
                f"Can't merge matrices with mixed shapes: {set(widths)}. Details: {list(zip(paths, widths))}"
            )
        rna = Counts.concatenate(rna_list)
        meta = None
        if metadata is not None:
            metadata_cols = [col for col in metadata.columns if not col.startswith("path_")]
            metadata = metadata[metadata_cols]
            cells_per_matrix = [counts.shape[0] for counts in rna_list]
            meta = metadata.loc[metadata.index.repeat(cells_per_matrix)].reset_index(drop=True)
            if id_col in metadata:
                rna.index = rna.index.str.slice(0, -1) + meta[id_col]
        if rna.index.duplicated().any():
            raise ValueError(
                "cell identifiers must be unique. Consider using metadata with `lane_id` column or specify a custom "
                "`id_col"
            )
        if meta is not None:
            meta.index = rna.cell_ids
        else:
            meta = rna.cell_ids
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            meta.to_csv(save_dir / "meta.tsv", sep="\t")
            # TODO: move create_rds val to config
            rna.save(save_dir / "rna.pickle", save_rds=True)
        return rna, meta

    @staticmethod
    def _merge_vdj(paths, metadata, save_dir):
        raise NotImplementedError()

    @staticmethod
    def _merge_surface(paths, metadata, save_dir):
        raise NotImplementedError()

    @staticmethod
    def _merge_antigen(paths, metadata, save_dir):
        raise NotImplementedError()

    @staticmethod
    def _merge_cnv(paths, metadata, save_dir):
        raise NotImplementedError()

    @staticmethod
    def _merge_atac(paths, metadata, save_dir):
        raise NotImplementedError()

    @staticmethod
    def _merge_spatial(paths, metadata, save_dir):
        raise NotImplementedError()

    @staticmethod
    def _merge_crispr(paths, metadata, save_dir):
        raise NotImplementedError()
