from dataforest.hooks import dataprocess


@dataprocess(requires="cluster")
def markers(branch: "CellBranch"):
    process_name = "markers"
    input_metadata_path = ProcessMethodsSC._get_temp_meta_path(branch, process_name)
    meta = branch.READER_METHODS.tsv(input_metadata_path, header=0)
    cluster_counts = meta["cluster_id"].value_counts()
    deficient_clusters = cluster_counts[cluster_counts < 2].index.tolist()
    if deficient_clusters:
        raise ValueError(
            f"Deficient clusters {deficient_clusters} with fewer than 3 "
            f"cells. Adjust clustering parameters.\n{cluster_counts}"
        )

    input_rds_path = branch["cluster"].path_map["cluster_r"]
    output_markers_path = branch["markers"].path
    logfc_thresh = branch.spec[process_name]["logfc_thresh"]
    r_functions_filepath = branch.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
    arg_list = [
        input_metadata_path,
        input_rds_path,
        output_markers_path,
        logfc_thresh,
        r_functions_filepath,
    ]
    r_find_markers = branch.schema.R_FILEPATHS["FIND_CLUSTER_MARKERS_SCRIPT"]
    ProcessMethodsSC._run_r_script(branch, r_find_markers, arg_list, "markers")


@dataprocess(requires="normalize", comparative=True)
def diffexp_bulk(branch: "CellBranch"):
    process_name = "diffexp_bulk"
    input_metadata_path = branch.get_temp_meta_path(branch, process_name)
    input_rds_path = branch["normalize"].path_map["matrix_r"]
    output_diffexp_path = branch[process_name].path_map["diffexp_bulk_result"]
    groups = branch[process_name].branch.meta["partition_code"].unique().astype("O")
    if len(groups) != 2:
        raise ValueError(f"Exactly two groups required for diffexp_bulk. Got: {groups}")
    test = branch.spec[process_name]["test"]
    ident1 = groups.min()
    ident2 = groups.max()
    groupby = "partition_code"
    logfc_thresh = branch.spec[process_name]["logfc_thresh"]
    r_functions_filepath = branch.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
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
    r_diff_exp_bulk_filepath = branch.schema.R_FILEPATHS["DIFF_EXP_BULK_SCRIPT"]
    ProcessMethodsSC._run_r_script(branch, r_diff_exp_bulk_filepath, arg_list, process_name)


@dataprocess(requires="cluster", comparative=True)
def diffexp(branch: "CellBranch"):
    # TODO: refactor both diffexp versions into `_get_diffexp_args`
    process_name = "diffexp"
    input_metadata_path = branch.get_temp_meta_path(branch, process_name)
    input_rds_path = branch["cluster"].path_map["cluster_r"]
    output_diffexp_path = branch[process_name].path_map["diffexp_result"]
    groups = branch[process_name].branch.meta["partition_code"].unique().astype("O")
    if len(groups) != 2:
        raise ValueError(f"Exactly two groups required for diffexp. Got: {groups}")
    test = branch.spec[process_name]["test"]
    ident1 = groups.min()
    ident2 = groups.max()
    groupby = "partition_code"
    logfc_thresh = branch.spec[process_name]["logfc_thresh"]
    r_functions_filepath = branch.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
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
    r_diff_exp_filepath = branch.schema.R_FILEPATHS["DIFF_EXP_CLUSTER_SCRIPT"]
    ProcessMethodsSC._run_r_script(branch, r_diff_exp_filepath, arg_list, process_name)
