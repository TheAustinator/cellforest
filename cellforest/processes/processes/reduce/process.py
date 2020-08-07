from pathlib import Path

from dataforest.hooks import dataprocess

from cellforest.processes.processes.reduce.umap import run_umap
from cellforest.utils.r.run_r_script import run_process_r_script
from cellforest.processes import R_FUNCTIONS_FILEPATH

R_PCA_SCRIPT = Path(__file__).parent / "pca.R"


@dataprocess(requires="normalize")
def reduce(branch: "CellBranch", run_name: str):
    input_metadata_path = branch.get_temp_meta_path(run_name)
    run_spec = branch.spec[run_name]
    params = run_spec.params
    process_run = branch[run_name]
    input_rds_path = process_run.path_map_prior["rna_r"]
    output_embeddings_path = process_run.path_map["pca_embeddings"]
    output_loadings_path = process_run.path_map["pca_loadings"]
    npcs = params["pca_npcs"]
    r_functions_filepath = R_FUNCTIONS_FILEPATH
    arg_list = [
        input_metadata_path,
        input_rds_path,
        output_embeddings_path,
        output_loadings_path,
        npcs,
        r_functions_filepath,
    ]
    run_process_r_script(branch, R_PCA_SCRIPT, arg_list, run_name)
    meta = run_umap(
        output_embeddings_path,
        n_neighbors=params["umap_n_neighbors"],
        min_dist=params["umap_min_dist"],
        n_components=params["umap_n_components"],
        metric=params["umap_metric"],
    )
    meta.index = branch.meta.index
    output_meta_path = process_run.path_map["meta"]
    meta.to_csv(output_meta_path, sep="\t")
    # branch.w["umap_embeddings"](meta, index=True, header=True)
