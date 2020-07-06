from dataforest.hooks import dataprocess


@dataprocess(requires="cluster")
def markers(forest: "CellForest"):
    process_name = "markers"
    input_metadata_path = ProcessMethodsSC._get_temp_meta_path(forest, process_name)
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


@dataprocess(requires="normalize", comparative=True)
def diffexp_bulk(forest: "CellForest"):
    process_name = "diffexp_bulk"
    input_metadata_path = forest.get_temp_meta_path(forest, process_name)
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


@dataprocess(requires="cluster", comparative=True)
def diffexp(forest: "CellForest"):
    # TODO: refactor both diffexp versions into `_get_diffexp_args`
    process_name = "diffexp"
    input_metadata_path = forest.get_temp_meta_path(forest, process_name)
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
