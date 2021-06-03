from anndata import AnnData


def prep(ad: AnnData) -> AnnData:
    """
    For compatibility with anndata2ri re-instantiates and converts all object columns to cat
    """
    ad = AnnData(ad.X, ad.obs, ad.var, ad.uns, ad.obsm, ad.varm, ad.layers, ad.raw, obsp=ad.obsp, varp=ad.varp)
    cols_obj = ad.obs.dtypes[ad.obs.dtypes == "object"].index.tolist()
    ad.obs[cols_obj] = ad.obs[cols_obj].astype("category")
    return ad
