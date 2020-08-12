from pathlib import Path

from dataforest.hooks import dataprocess

from cellforest.processes import R_FUNCTIONS_FILEPATH
from cellforest.utils.r.run_r_script import run_process_r_script

R_MARKERS_SCRIPT = Path(__file__).parent / "markers.R"
R_DIFFEXP_SCRIPT = Path(__file__).parent / "diffexp.R"


@dataprocess(requires="cluster")
def markers(branch: "CellBranch", run_name: str):
    input_metadata_path = branch.get_temp_meta_path(run_name)
    run_spec = branch.spec[run_name]
    params = run_spec.params
    process_run = branch[run_name]
    output_markers_path = process_run.path
    root_dir = str(branch.root)
    spec_str = branch.spec.shell_str
    cluster_counts = branch.meta["cluster_id"].value_counts()
    deficient_clusters = cluster_counts[cluster_counts < 2].index.tolist()
    if deficient_clusters:
        raise ValueError(
            f"Deficient clusters {deficient_clusters} with fewer than 3 "
            f"cells. Adjust clustering parameters.\n{cluster_counts}"
        )
    arg_list = [
        input_metadata_path,
        output_markers_path,
        root_dir,
        spec_str,
        params["logfc_thresh"],
        params["test"],
        R_FUNCTIONS_FILEPATH,
    ]
    run_process_r_script(branch, R_MARKERS_SCRIPT, arg_list, run_name)


@dataprocess(requires="cluster", comparative=True)
def diffexp(branch: "CellBranch", run_name: str):
    # TODO: refactor both diffexp versions into `_get_diffexp_args`
    input_metadata_path = branch.get_temp_meta_path(run_name)
    run_spec = branch.spec[run_name]
    params = run_spec.params
    process_run = branch[run_name]
    output_diffexp_path = process_run.path_map["diffexp"]
    root_dir = str(branch.root)
    spec_str = branch.spec.shell_str
    groups = branch[run_name].branch.meta["partition_code"].unique().astype("O")
    if len(groups) != 2:
        raise ValueError(f"Exactly two groups required for diffexp. Got: {groups}")
    ident1 = groups.min()
    ident2 = groups.max()
    groupby = "partition_code"
    arg_list = [
        input_metadata_path,
        output_diffexp_path,
        root_dir,
        spec_str,
        params["test"],
        params["logfc_thresh"],
        ident1,
        ident2,
        groupby,
        R_FUNCTIONS_FILEPATH,
    ]
    run_process_r_script(branch, R_DIFFEXP_SCRIPT, arg_list, run_name)


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
