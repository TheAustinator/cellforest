from dataforest.hooks import dataprocess


@dataprocess(requires="normalize", comparative=True, temp_meta=False)
def gsea_bulk(forest: "CellBranch", **kwargs):
    # from scgsea.GSEA import GSEA

    gsea = GSEA(forest.at("gsea"), "gsea")
    gsea.run(**kwargs)
    return gsea


@dataprocess(requires="cluster", comparative=True, temp_meta=False)
def gsea(forest: "CellBranch"):
    from gsea.GSEA import GSEA

    gsea = GSEA(forest["gsea"].branch, "gsea")
    gsea.run()
