from dataforest.hooks import dataprocess


@dataprocess(requires="normalize", comparative=True, temp_meta=False)
def gsea_bulk(forest: "CellForest", **kwargs):
    # from scgsea.GSEA import GSEA

    gsea = GSEA(forest.at("gsea"), "gsea")
    gsea.run(**kwargs)
    return gsea


@dataprocess(requires="cluster", comparative=True, temp_meta=False)
def gsea(forest: "CellForest"):
    from gsea.GSEA import GSEA

    gsea = GSEA(forest["gsea"].forest, "gsea")
    gsea.run()
