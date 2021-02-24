import scanpy as sc


def add_var_pos(var):
    """Add chromosomal center position AnnData var"""
    pos = sc.queries.biomart_annotations(
        "hsapiens", ["ensembl_gene_id", "start_position", "end_position", "chromosome_name"]
    )
    pos.set_index("ensembl_gene_id", inplace=True)
    pos["center_position"] = pos[["start_position", "end_position"]].mean(axis=1)
    pos = pos["center_position"]
    var = var.merge(pos, left_on="ensgs", right_index=True, how="left")
    return var
