from dataforest.hooks import hook

from cellforest.utils.r.Convert import Convert


@hook(attrs=["matrix_layer"])
def hook_unify_matrix_node(dp):
    """
    If node has counts matrix output, ensure that all desired formats are
    present (e.g. Counts pickle, Seurat rds, anndata, cellranger). If any are
    missing, run the necessary conversions.
    Note: only Counts pickle and Seurat rds currently implemented
    """
    # TODO: currently hardcoded to pickle and rds, but fix this later
    if dp.matrix_layer:
        pickle_path = dp.branch[dp.name].path_map["rna"]
        rds_path = dp.branch[dp.name].path_map["rna_r"]
        if dp.output == "pickle":
            Convert.pickle_to_rds_dir(pickle_path.parent)
        elif dp.output == "rds":
            Convert.rds_to_pickle_dir(rds_path.parent)
        else:
            raise NotImplementedError(f"Only rds and pickle outputs currently supported")
