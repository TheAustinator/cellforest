from pathlib import Path

from dataforest.hooks import dataprocess

# TODO: what to do about core/utility methods? core module? move to utils?
from cellforest.utils.r.run_r_script import run_process_r_script
from cellforest.processes import R_FUNCTIONS_FILEPATH

R_SCTRANSFORM_SCRIPT = Path(__file__).parent / "sctransform.R"
R_SEURAT_DEFAULT_NORM_SCRIPT = Path(__file__).parent / "seurat_default_normalize.R"


@dataprocess(matrix_layer=True, output="rds")
def normalize(branch: "CellBranch", run_name: str):
    """
    Performs:
        - cell filtering by `min_genes`, `max_genes`, and `perc_mito_cutoff`
        - gene filtering by `min_cells` expressing
        - normalization via either:
            - seurat default
            - sctransform
    Params:
        min_genes (int):
        max_genes (int):
        min_cells (int):
        perc_mito_cutoff (int, float):
        method (str): from {"seurat_default", "sctransform"}
        nfeatures (int): (seurat_default only)
    """
    input_metadata_path = branch.get_temp_meta_path(run_name)
    # TODO: add a root filepaths lookup
    run_spec = branch.spec[run_name]
    params = run_spec.params
    process_run = branch[run_name]
    input_rds_path = process_run.path_map_prior["rna_r"]
    output_rds_path = process_run.path_map["rna_r"]
    min_genes = params["min_genes"]
    max_genes = params["max_genes"]
    min_cells = params["min_cells"]
    perc_mito_cutoff = params["perc_mito_cutoff"]
    method = params["method"]
    arg_list = [
        input_metadata_path,
        input_rds_path,
        output_rds_path,
        min_genes,
        max_genes,
        min_cells,
        perc_mito_cutoff,
        R_FUNCTIONS_FILEPATH,
    ]
    if method == "sctransform":
        output_corrected_umi_path = process_run.path_map["corrected_umi"]
        output_pearson_residual_path = process_run.path_map["pearson_residual"]
        arg_list += [output_corrected_umi_path, output_pearson_residual_path]
        r_normalize_script = R_SCTRANSFORM_SCRIPT
    elif method == "seurat_default":
        verbose = True
        verbose = str(verbose).upper()
        nfeatures = params["nfeatures"]
        arg_list += [verbose, nfeatures]
        r_normalize_script = R_SEURAT_DEFAULT_NORM_SCRIPT
    else:
        raise ValueError(f"Invalid normalization method: {method}. Use 'sctransform' or 'seurat_default'")
    run_process_r_script(branch, r_normalize_script, arg_list, run_name)
