from typing import Union, TYPE_CHECKING

from dataforest.templates.ProcessMethods import ProcessMethods
from pathlib import Path
import pandas as pd

from cellforest.utils.r.seurat_rds_to_pickle import seurat_rds_to_sparse_pickle
from cellforest.utils.r.shell_command import process_shell_command
from cellforest.templates.dataprocess_sc import dataprocess_sc

if TYPE_CHECKING:
    from cellforest.templates.CellForest import CellForest


class ProcessMethodsSC(ProcessMethods):
    @staticmethod
    @dataprocess_sc(requires="root")
    def normalize(forest: "CellForest"):
        process_name = "normalize"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(forest, process_name)
        input_dir = forest.root_dir
        # TODO: current position
        output_rds_path = forest[process_name].path_map["matrix_r"]
        min_genes = forest.spec[process_name]["min_genes"]
        max_genes = forest.spec[process_name]["max_genes"]
        min_cells = forest.spec[process_name]["min_cells"]
        perc_mito_cutoff = forest.spec[process_name]["perc_mito_cutoff"]
        r_functions_filepath = forest.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        method = forest.spec[process_name]["method"]
        arg_list = [
            input_metadata_path,
            input_dir,
            output_rds_path,
            min_genes,
            max_genes,
            min_cells,
            perc_mito_cutoff,
            r_functions_filepath,
        ]
        if method == "sctransform":
            output_corrected_umi_path = forest[process_name].path_map["corrected_umi"]
            output_pearson_residual_path = forest[process_name].path_map["pearson_residual"]
            arg_list += [output_corrected_umi_path, output_pearson_residual_path]
            r_normalize_script = str(forest.schema.R_FILEPATHS["SCTRANSFORM_SCRIPT"])
        elif method == "seurat_default":
            verbose = True
            verbose = str(verbose).upper()
            nfeatures = forest.spec[process_name]["nfeatures"]
            arg_list += [verbose, nfeatures]
            r_normalize_script = str(forest.schema.R_FILEPATHS["SEURAT_DEFAULT_NORMALIZE_SCRIPT"])
        else:
            raise ValueError(f"Invalid normalization method: {method}. Use 'sctransform' or 'seurat_default'")
        ProcessMethodsSC._run_r_script(forest, r_normalize_script, arg_list, process_name)
        seurat_rds_to_sparse_pickle(forest.paths["combine"], output_rds_path, forest.paths[process_name])

    @staticmethod
    @dataprocess_sc(requires="normalize", comparative=True, temp_meta=False)
    def gsea_bulk(forest: "CellForest", **kwargs):
        # from scgsea.GSEA import GSEA

        gsea = GSEA(forest.at("gsea_bulk"), "gsea_bulk")
        gsea.run(**kwargs)
        return gsea

    @staticmethod
    @dataprocess_sc(requires="normalize")
    def dim_reduce(forest: "CellForest"):
        process_name = "dim_reduce"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(forest, process_name)
        input_rds_path = forest["normalize"].path_map["matrix_r"]
        print(input_rds_path)
        output_rds_path = forest[process_name].path_map["dimred_r"]
        output_embeddings_path = forest[process_name].path_map["pca_embeddings"]
        output_loadings_path = forest[process_name].path_map["pca_loadings"]
        npcs = forest.spec[process_name]["pca_npcs"]
        r_functions_filepath = forest.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        arg_list = [
            input_metadata_path,
            input_rds_path,
            output_rds_path,
            output_embeddings_path,
            output_loadings_path,
            npcs,
            r_functions_filepath,
        ]
        r_pca_filepath = str(forest.schema.R_FILEPATHS["PCA_SCRIPT"])
        working_dir = str(forest.paths[process_name])
        arg_list = list(map(str, arg_list))
        command_string = f"Rscript {r_pca_filepath} {' '.join(arg_list)}"
        process_shell_command(command_string=command_string, working_dir=working_dir, process_name="pca")
        umap_df = ProcessMethodsSC._run_umap(
            output_embeddings_path,
            n_neighbors=forest.spec[process_name]["umap_n_neighbors"],
            min_dist=forest.spec[process_name]["umap_min_dist"],
            n_components=forest.spec[process_name]["umap_n_components"],
            metric=forest.spec[process_name]["umap_metric"],
        )
        umap_df.index = forest.f["normalize"]["cell_ids"][0]
        output_umap_embeddings_path = forest[process_name].path_map["umap_embeddings"]
        forest.write_umap_embeddings(output_umap_embeddings_path, umap_df, index=True, header=True)

    @staticmethod
    @dataprocess_sc(requires="dim_reduce")
    def cluster(forest: "CellForest"):
        process_name = "cluster"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(forest, process_name)
        input_rds_path = forest["dim_reduce"].path_map["dimred_r"]
        output_rds_path = forest[process_name].path_map["cluster_r"]
        output_clusters_path = forest[process_name].path_map["clusters"]
        num_pcs = forest.spec[process_name]["num_pcs"]
        res = forest.spec[process_name]["res"]
        eps = forest.spec[process_name]["eps"]
        r_functions_filepath = forest.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
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
        r_clusters_filepath = forest.schema.R_FILEPATHS["FIND_CLUSTERS_SCRIPT"]
        ProcessMethodsSC._run_r_script(forest, r_clusters_filepath, arg_list, process_name)

    @staticmethod
    @dataprocess_sc(requires="normalize", comparative=True)
    def diffexp_bulk(forest: "CellForest"):
        process_name = "diffexp_bulk"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(forest, process_name)
        input_rds_path = forest["normalize"].path_map["matrix_r"]
        output_diffexp_path = forest[process_name].path_map["diffexp_bulk_result"]
        groups = forest[process_name].forest.meta["partition_code"].unique().astype("O")
        if len(groups) != 2:
            raise ValueError(f"Exactly two groups required for diffexp_bulk. Got: {groups}")
        test = forest.spec[process_name]["test"]
        ident1 = groups.min()
        ident2 = groups.max()
        groupby = "partition_code"
        logfc_thresh = forest.spec[process_name]["logfc_thresh"]
        r_functions_filepath = forest.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
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
        r_diff_exp_bulk_filepath = forest.schema.R_FILEPATHS["DIFF_EXP_BULK_SCRIPT"]
        ProcessMethodsSC._run_r_script(forest, r_diff_exp_bulk_filepath, arg_list, process_name)

    @staticmethod
    @dataprocess_sc(requires="cluster", comparative=True)
    def diffexp(forest: "CellForest"):
        # TODO: refactor both diffexp versions into `_get_diffexp_args`
        process_name = "diffexp"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(forest, process_name)
        input_rds_path = forest["cluster"].path_map["cluster_r"]
        output_diffexp_path = forest[process_name].path_map["diffexp_result"]
        groups = forest[process_name].forest.meta["partition_code"].unique().astype("O")
        if len(groups) != 2:
            raise ValueError(f"Exactly two groups required for diffexp. Got: {groups}")
        test = forest.spec[process_name]["test"]
        ident1 = groups.min()
        ident2 = groups.max()
        groupby = "partition_code"
        logfc_thresh = forest.spec[process_name]["logfc_thresh"]
        r_functions_filepath = forest.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
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
        r_diff_exp_filepath = forest.schema.R_FILEPATHS["DIFF_EXP_CLUSTER_SCRIPT"]
        ProcessMethodsSC._run_r_script(forest, r_diff_exp_filepath, arg_list, process_name)

    @staticmethod
    @dataprocess_sc(requires="cluster", comparative=True, temp_meta=False)
    def gsea(forest: "CellForest"):
        from gsea.GSEA import GSEA

        gsea = GSEA(forest["gsea"].forest, "gsea")
        gsea.run()

    @staticmethod
    @dataprocess_sc(requires="cluster")
    def markers(forest: "CellForest"):
        process_name = "markers"
        input_metadata_path = ProcessMethodsSC._get_temp_metadata_path(forest, process_name)
        meta = forest.READER_METHODS.tsv(input_metadata_path, header=0)
        cluster_counts = meta["cluster_id"].value_counts()
        deficient_clusters = cluster_counts[cluster_counts < 2].index.tolist()
        if deficient_clusters:
            raise ValueError(
                f"Deficient clusters {deficient_clusters} with fewer than 3 "
                f"cells. Adjust clustering parameters.\n{cluster_counts}"
            )

        input_rds_path = forest["cluster"].path_map["cluster_r"]
        output_markers_path = forest["markers"].path
        logfc_thresh = forest.spec[process_name]["logfc_thresh"]
        r_functions_filepath = forest.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
        arg_list = [
            input_metadata_path,
            input_rds_path,
            output_markers_path,
            logfc_thresh,
            r_functions_filepath,
        ]
        r_find_markers = forest.schema.R_FILEPATHS["FIND_CLUSTER_MARKERS_SCRIPT"]
        ProcessMethodsSC._run_r_script(forest, r_find_markers, arg_list, "markers")

    @staticmethod
    def _run_r_script(forest: "CellForest", r_script_filepath: str, arg_list: list, process_name: str):
        command_string = f"Rscript {r_script_filepath} {' '.join(map(str, arg_list))}"
        working_dir = str(forest[process_name].path)
        process_shell_command(
            command_string=command_string, working_dir=working_dir, process_name=process_name,
        )

    @staticmethod
    def _get_temp_metadata_path(forest: "CellForest", process_name: str):
        return forest[process_name].path / forest.schema.TEMP_METADATA_FILENAME

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
            n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed, n_components=n_components, metric=metric,
        )
        umap_matrix = umap_handle.fit(pca_df).embedding_
        umap_df = pd.DataFrame(umap_matrix, columns=[f"UMAP_{idx + 1}" for idx in range(umap_matrix.shape[1])],)
        return umap_df
