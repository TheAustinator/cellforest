from dataforest.hooks import dataprocess


@dataprocess(requires="reduce")
def cluster(branch: "CellBranch"):
    process_name = "cluster"
    input_metadata_path = branch.get_temp_meta_path(branch, process_name)
    input_rds_path = branch["reduce"].path_map["dimred_r"]
    output_rds_path = branch[process_name].path_map["cluster_r"]
    output_clusters_path = branch[process_name].path_map["clusters"]
    num_pcs = branch.spec[process_name]["num_pcs"]
    res = branch.spec[process_name]["res"]
    eps = branch.spec[process_name]["eps"]
    r_functions_filepath = branch.schema.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
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
    r_clusters_filepath = branch.schema.R_FILEPATHS["FIND_CLUSTERS_SCRIPT"]
    ProcessMethodsSC._run_r_script(branch, r_clusters_filepath, arg_list, process_name)
