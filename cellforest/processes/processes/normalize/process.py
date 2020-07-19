from dataforest.hooks import dataprocess

# TODO: what to do about core/utility methods? core module? move to utils?
from cellforest.utils.r.run_r_script import run_process_r_script


@dataprocess(matrix_layer=True)
def normalize(forest: "CellForest", run_name):
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
    input_metadata_path = forest.get_temp_meta_path(run_name)
    # TODO: add a root filepaths lookup
    params = forest.spec[run_name].params
    input_rds_path = forest.root_dir / "rna.rds"
    output_rds_path = forest[run_name].path_map["rna_r"]
    min_genes = params["min_genes"]
    max_genes = params["max_genes"]
    min_cells = params["min_cells"]
    perc_mito_cutoff = params["perc_mito_cutoff"]
    r_functions_filepath = forest.schema.__class__.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
    method = params["method"]
    arg_list = [
        input_metadata_path,
        input_rds_path,
        output_rds_path,
        min_genes,
        max_genes,
        min_cells,
        perc_mito_cutoff,
        r_functions_filepath,
    ]
    if method == "sctransform":
        output_corrected_umi_path = forest[run_name].path_map["corrected_umi"]
        output_pearson_residual_path = forest[run_name].path_map["pearson_residual"]
        arg_list += [output_corrected_umi_path, output_pearson_residual_path]
        r_normalize_script = str(forest.schema.__class__.R_FILEPATHS["SCTRANSFORM_SCRIPT"])
    elif method == "seurat_default":
        verbose = True
        verbose = str(verbose).upper()
        nfeatures = params["nfeatures"]
        arg_list += [verbose, nfeatures]
        r_normalize_script = str(forest.schema.__class__.R_FILEPATHS["SEURAT_DEFAULT_NORMALIZE_SCRIPT"])
    else:
        raise ValueError(f"Invalid normalization method: {method}. Use 'sctransform' or 'seurat_default'")
    run_process_r_script(forest, r_normalize_script, arg_list, run_name)
