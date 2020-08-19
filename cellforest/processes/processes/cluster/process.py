from pathlib import Path

from dataforest.hooks import dataprocess

from cellforest.processes import R_FUNCTIONS_FILEPATH
from cellforest.utils.r.run_r_script import run_process_r_script

R_CLUSTER_SCRIPT = Path(__file__).parent / "cluster.R"


@dataprocess(requires="reduce")
def cluster(branch: "CellBranch", run_name: str):
    input_metadata_path = branch.get_temp_meta_path(run_name)
    # TODO: add a root filepaths lookup
    run_spec = branch.spec[run_name]
    params = run_spec.params
    process_run = branch[run_name]
    output_clusters_path = process_run.path_map["meta"]
    # TODO: may actually need way to pass spec through
    root_dir = str(branch.root)
    spec_str = branch.spec.shell_str
    arg_list = [
        input_metadata_path,
        output_clusters_path,
        root_dir,
        spec_str,
        params["num_pcs"],
        params["res"],
        params["eps"],
        R_FUNCTIONS_FILEPATH,
    ]
    run_process_r_script(branch, R_CLUSTER_SCRIPT, arg_list, run_name)
