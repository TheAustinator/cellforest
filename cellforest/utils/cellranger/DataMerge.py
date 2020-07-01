import os

from cellforest import Counts


class DataMerge:
    @staticmethod
    def merge_assay(paths, mode, metadata=None, save_dir=None):
        method = getattr(DataMerge, f"_merge_{mode}")
        return method(paths, metadata, save_dir)

    @staticmethod
    def _merge_rna(paths, metadata, save_dir):
        """"""
        rna_list = [Counts.from_cellranger(dir_) for dir_ in paths]
        meta = None
        if metadata is not None:
            metadata_cols = [col for col in metadata.columns if not col.startswith("path_")]
            metadata = metadata[metadata_cols]
            cells_per_matrix = [counts.shape[0] for counts in rna_list]
            meta = metadata.loc[metadata.index.repeat(cells_per_matrix)].reset_index(drop=True)
        rna = Counts.concatenate(rna_list)
        if meta is not None:
            meta.index = rna.cell_ids
        else:
            meta = rna.cell_ids
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            meta.to_csv(save_dir / "meta.tsv", sep="\t")
            # TODO: move create_rds val to config
            rna.save(save_dir / "rna.pickle", create_rds=True)
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
