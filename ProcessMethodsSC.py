from pathlib import Path
from typing import Union, TYPE_CHECKING

import pandas as pd

from dataforest.ProcessMethods import ProcessMethods
from cellforest.r_utils.seurat_rds_to_pickle import seurat_rds_to_sparse_pickle
from cellforest.r_utils.shell_command import shell_command
from cellforest.dataprocess_sc import dataprocess_sc

if TYPE_CHECKING:
    from cellforest.CellORM import CellForest


class ProcessMethodsSC(ProcessMethods):
    def __init__(self):
        self.normalize = self.normalize()

    @staticmethod
    def combine(root_path: Union[str, Path], sample_metadata_df: pd.DataFrame, **kwargs):
        from cellforest.preprocessing import combine_10x_outputs

        path = Path(root_path)
        if not path.exists():
            path.mkdir(parents=True, exist_ok=True)
        return combine_10x_outputs(root_path, sample_metadata_df, **kwargs)

    @staticmethod
    @dataprocess_sc(requires="combine")
    def normalize(orm: "CellForest"):
        process_name = "normalize"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(orm, process_name)
        input_tenx_directory_path = orm["combine"].path
        output_rds_path = orm[process_name].path_map["matrix_r"]
        min_genes = orm.spec[process_name]["min_genes"]
        max_genes = orm.spec[process_name]["max_genes"]
        min_cells = orm.spec[process_name]["min_cells"]
        perc_mito_cutoff = orm.spec[process_name]["perc_mito_cutoff"]
        r_functions_filepath = orm.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        method = orm.spec[process_name]["method"]
        arg_list = [
            input_metadata_path,
            input_tenx_directory_path,
            output_rds_path,
            min_genes,
            max_genes,
            min_cells,
            perc_mito_cutoff,
            r_functions_filepath,
        ]
        if method == "sctransform":
            output_corrected_umi_path = orm[process_name].path_map["corrected_umi"]
            output_pearson_residual_path = orm[process_name].path_map["pearson_residual"]
            arg_list += [output_corrected_umi_path, output_pearson_residual_path]
            r_normalize_script = str(orm.schema.R_FILEPATHS["SCTRANSFORM_SCRIPT"])
        elif method == "seurat_default":
            verbose = True
            verbose = str(verbose).upper()
            nfeatures = orm.spec[process_name]["nfeatures"]
            arg_list += [verbose, nfeatures]
            r_normalize_script = str(orm.schema.R_FILEPATHS["SEURAT_DEFAULT_NORMALIZE_SCRIPT"])
        else:
            raise ValueError(f"Invalid normalization method: {method}. Use 'sctransform' or 'seurat_default'")
        ProcessMethodsSC._run_r_script(orm, r_normalize_script, arg_list, process_name)
        seurat_rds_to_sparse_pickle(orm.paths["combine"], output_rds_path, orm.paths[process_name])

    @staticmethod
    @dataprocess_sc(requires="normalize", comparative=True, temp_meta=False)
    def gsea_bulk(orm: "CellForest", **kwargs):
        from gsea.GSEA import GSEA

        gsea = GSEA(orm.at("gsea_bulk"), "gsea_bulk")
        gsea.run(**kwargs)
        return gsea

    @staticmethod
    @dataprocess_sc(requires="normalize")
    def dim_reduce(orm: "CellForest"):
        process_name = "dim_reduce"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(orm, process_name)
        input_rds_path = orm["normalize"].path_map["matrix_r"]
        print(input_rds_path)
        output_rds_path = orm[process_name].path_map["dimred_r"]
        output_embeddings_path = orm[process_name].path_map["pca_embeddings"]
        output_loadings_path = orm[process_name].path_map["pca_loadings"]
        npcs = orm.spec[process_name]["pca_npcs"]
        r_functions_filepath = orm.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        arg_list = [
            input_metadata_path,
            input_rds_path,
            output_rds_path,
            output_embeddings_path,
            output_loadings_path,
            npcs,
            r_functions_filepath,
        ]
        r_pca_filepath = str(orm.schema.R_FILEPATHS["PCA_SCRIPT"])
        working_dir = str(orm.paths[process_name])
        arg_list = list(map(str, arg_list))
        command_string = f"Rscript {r_pca_filepath} {' '.join(arg_list)}"
        shell_command(command_string=command_string, working_dir=working_dir, process_name="pca")
        umap_df = ProcessMethodsSC._run_umap(
            output_embeddings_path,
            n_neighbors=orm.spec[process_name]["umap_n_neighbors"],
            min_dist=orm.spec[process_name]["umap_min_dist"],
            n_components=orm.spec[process_name]["umap_n_components"],
            metric=orm.spec[process_name]["umap_metric"],
        )
        umap_df.index = orm.f_cell_ids[0]
        output_umap_embeddings_path = orm[process_name].path_map["umap_embeddings"]
        orm.write_umap_embeddings(output_umap_embeddings_path, umap_df, index=True, header=True)

    @staticmethod
    @dataprocess_sc(requires="dim_reduce")
    def cluster(orm: "CellForest"):
        process_name = "cluster"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(orm, process_name)
        input_rds_path = orm["dim_reduce"].path_map["dimred_r"]
        output_rds_path = orm[process_name].path_map["cluster_r"]
        output_clusters_path = orm[process_name].path_map["clusters"]
        num_pcs = orm.spec[process_name]["num_pcs"]
        res = orm.spec[process_name]["res"]
        eps = orm.spec[process_name]["eps"]
        r_functions_filepath = orm.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        arg_list = [
            input_metadata_path,
            input_rds_path,
            output_rds_path,
            output_clusters_path,
            num_pcs,
            res,
            eps,
            r_functions_filepath,
        ]
        r_clusters_filepath = orm.schema.R_FILEPATHS["FIND_CLUSTERS_SCRIPT"]
        ProcessMethodsSC._run_r_script(orm, r_clusters_filepath, arg_list, process_name)

    @staticmethod
    @dataprocess_sc(requires="normalize", comparative=True)
    def diffexp_bulk(orm: "CellForest"):
        process_name = "diffexp_bulk"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(orm, process_name)
        input_rds_path = orm["normalize"].path_map["matrix_r"]
        output_diffexp_path = orm[process_name].path_map["diffexp_bulk_result"]
        groups = orm[process_name].orm.meta["partition_code"].unique().astype("O")
        if len(groups) != 2:
            raise ValueError(f"Exactly two groups required for diffexp_bulk. Got: {groups}")
        test = orm.spec[process_name]["test"]
        ident1 = groups.min()
        ident2 = groups.max()
        groupby = "partition_code"
        logfc_thresh = orm.spec[process_name]["logfc_thresh"]
        r_functions_filepath = orm.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        arg_list = [
            input_metadata_path,
            input_rds_path,
            output_diffexp_path,
            test,
            ident1,
            ident2,
            groupby,
            logfc_thresh,
            r_functions_filepath,
        ]
        r_diff_exp_bulk_filepath = orm.schema.R_FILEPATHS["DIFF_EXP_BULK_SCRIPT"]
        ProcessMethodsSC._run_r_script(orm, r_diff_exp_bulk_filepath, arg_list, process_name)

    @staticmethod
    @dataprocess_sc(requires="cluster", comparative=True)
    def diffexp(orm: "CellForest"):
        # TODO: refactor both diffexp versions into `_get_diffexp_args`
        process_name = "diffexp"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(orm, process_name)
        input_rds_path = orm["cluster"].path_map["cluster_r"]
        output_diffexp_path = orm[process_name].path_map["diffexp_result"]
        groups = orm[process_name].orm.meta["partition_code"].unique().astype("O")
        if len(groups) != 2:
            raise ValueError(f"Exactly two groups required for diffexp. Got: {groups}")
        test = orm.spec[process_name]["test"]
        ident1 = groups.min()
        ident2 = groups.max()
        groupby = "partition_code"
        logfc_thresh = orm.spec[process_name]["logfc_thresh"]
        r_functions_filepath = orm.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        arg_list = [
            input_metadata_path,
            input_rds_path,
            output_diffexp_path,
            test,
            ident1,
            ident2,
            groupby,
            logfc_thresh,
            r_functions_filepath,
        ]
        r_diff_exp_filepath = orm.schema.R_FILEPATHS["DIFF_EXP_CLUSTER_SCRIPT"]
        ProcessMethodsSC._run_r_script(orm, r_diff_exp_filepath, arg_list, process_name)

    @staticmethod
    @dataprocess_sc(requires="cluster", comparative=True, temp_meta=False)
    def gsea(orm: "CellForest"):
        from gsea.GSEA import GSEA

        gsea = GSEA(orm["gsea"].orm, "gsea")
        gsea.run()

    @staticmethod
    @dataprocess_sc(requires="cluster")
    def markers(orm: "CellForest"):
        process_name = "markers"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(orm, process_name)
        meta = orm.READER_METHODS.tsv(input_metadata_path, header=0)
        cluster_counts = meta["cluster_id"].value_counts()
        deficient_clusters = cluster_counts[cluster_counts < 2].index.tolist()
        if deficient_clusters:
            raise ValueError(
                f"Deficient clusters {deficient_clusters} with fewer than 3 "
                f"cells. Adjust clustering parameters.\n{cluster_counts}"
            )

        input_rds_path = orm["cluster"].path_map["cluster_r"]
        output_markers_path = orm["markers"].path
        logfc_thresh = orm.spec[process_name]["logfc_thresh"]
        r_functions_filepath = orm.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        arg_list = [input_metadata_path, input_rds_path, output_markers_path, logfc_thresh, r_functions_filepath]
        r_find_markers = orm.schema.R_FILEPATHS["FIND_CLUSTER_MARKERS_SCRIPT"]
        ProcessMethodsSC._run_r_script(orm, r_find_markers, arg_list, "markers")

    @staticmethod
    def _run_r_script(orm: "CellForest", r_script_filepath: str, arg_list: list, process_name: str):
        command_string = f"Rscript {r_script_filepath} {' '.join(map(str, arg_list))}"
        working_dir = str(orm[process_name].path)
        shell_command(command_string=command_string, working_dir=working_dir, process_name=process_name)

    @staticmethod
    def _get_temp_metadata_path(orm: "CellForest", process_name: str):
        return orm[process_name].path / orm.schema.TEMP_METADATA_FILENAME

    @staticmethod
    def _run_umap(
        input_file_path: str,
        n_neighbors: int = 10,
        min_dist: float = 0.5,
        n_components: int = 2,
        metric: str = "euclidean",
        seed: int = 42,
    ):
        import umap
        import warnings
        from numba.errors import NumbaPerformanceWarning

        warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)

        pca_df = pd.read_csv(input_file_path, sep="\t")
        umap_handle = umap.UMAP(
            n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed, n_components=n_components, metric=metric
        )
        umap_matrix = umap_handle.fit(pca_df).embedding_
        umap_df = pd.DataFrame(umap_matrix, columns=[f"UMAP_{idx + 1}" for idx in range(umap_matrix.shape[1])])
        return umap_df
