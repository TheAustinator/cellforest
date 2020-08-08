from pathlib import Path

from dataforest.hooks import dataprocess

from cellforest.processes import R_FUNCTIONS_FILEPATH
from cellforest.utils.r.run_r_script import run_process_r_script

R_SCTRANSFORM_SCRIPT = Path(__file__).parent / "cluster.R"


@dataprocess(requires="reduce")
def cluster(branch: "CellBranch", run_name: str):
    input_metadata_path = branch.get_temp_meta_path(run_name)
    # TODO: add a root filepaths lookup
    run_spec = branch.spec[run_name]
    params = run_spec.params
    process_run = branch[run_name]

    input_rds_path = process_run.path_map["matrix_r"]
    output_rds_path = process_run.path_map["cluster_r"]
    output_clusters_path = process_run.path_map["clusters"]
    num_pcs = params["num_pcs"]
    res = params["res"]
    eps = params["eps"]
    arg_list = [
        input_metadata_path,
        input_rds_path,
        output_rds_path,
        output_clusters_path,
        num_pcs,
        res,
        eps,
        R_FUNCTIONS_FILEPATH,
    ]
    r_clusters_filepath = branch.schema.R_FILEPATHS["FIND_CLUSTERS_SCRIPT"]
    run_process_r_script(branch, r_clusters_filepath, arg_list, run_name)
