from dataforest.hooks import dataprocess


# TODO: strange bug where normalize attempts to run this rather than normalize, so commented out
@dataprocess(requires="normalize")
def reduce(forest: "CellForest"):
    process_name = "reduce"
    input_metadata_path = forest.get_temp_metadata_path(forest, process_name)
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


# def _run_umap(
#     input_file_path: str,
#     n_neighbors: int = 10,
#     min_dist: float = 0.5,
#     n_components: int = 2,
#     metric: str = "euclidean",
#     seed: int = 42,
# ):
#     import umap
#     import warnings
#     from numba.errors import NumbaPerformanceWarning
#
#     warnings.filterwarnings("ignore", category=NumbaPerformanceWarning)
#
#     pca_df = pd.read_csv(input_file_path, sep="\t")
#     umap_handle = umap.UMAP(
#         n_neighbors=n_neighbors, min_dist=min_dist, random_state=seed, n_components=n_components, metric=metric,
#     )
#     umap_matrix = umap_handle.fit(pca_df).embedding_
#     umap_df = pd.DataFrame(umap_matrix, columns=[f"UMAP_{idx + 1}" for idx in range(umap_matrix.shape[1])],)
#     return umap_df
