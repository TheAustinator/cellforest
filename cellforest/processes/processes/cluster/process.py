from dataforest.hooks import dataprocess


@dataprocess(requires="dim_reduce")
def cluster(forest: "CellForest"):
    process_name = "cluster"
    input_metadata_path = forest.get_temp_metadata_path(forest, process_name)
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
