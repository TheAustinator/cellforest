from dataforest.hooks import dataprocess

# TODO: what to do about core/utility methods? core module? move to utils?
from cellforest.utils.r.run_r_script import run_process_r_script


@dataprocess(requires="root", matrix_layer=True)
def normalize(forest: "CellForest"):
    process_name = "normalize"
    input_metadata_path = forest.get_temp_metadata_path(process_name)
    # TODO: add a root filepaths lookup
    input_rds_path = forest.root_dir / "rna.rds"
    output_rds_path = forest[process_name].path_map["rna_r"]
    min_genes = forest.spec[process_name]["min_genes"]
    max_genes = forest.spec[process_name]["max_genes"]
    min_cells = forest.spec[process_name]["min_cells"]
    perc_mito_cutoff = forest.spec[process_name]["perc_mito_cutoff"]
    r_functions_filepath = forest.schema.__class__.R_FILEPATHS["FUNCTIONS_FILE_PATH"]
    method = forest.spec[process_name]["method"]
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
        output_corrected_umi_path = forest[process_name].path_map["corrected_umi"]
        output_pearson_residual_path = forest[process_name].path_map["pearson_residual"]
        arg_list += [output_corrected_umi_path, output_pearson_residual_path]
        r_normalize_script = str(forest.schema.__class__.R_FILEPATHS["SCTRANSFORM_SCRIPT"])
    elif method == "seurat_default":
        verbose = True
        verbose = str(verbose).upper()
        nfeatures = forest.spec[process_name]["nfeatures"]
        arg_list += [verbose, nfeatures]
        r_normalize_script = str(forest.schema.__class__.R_FILEPATHS["SEURAT_DEFAULT_NORMALIZE_SCRIPT"])
    else:
        raise ValueError(f"Invalid normalization method: {method}. Use 'sctransform' or 'seurat_default'")
    run_process_r_script(forest, r_normalize_script, arg_list, process_name)
