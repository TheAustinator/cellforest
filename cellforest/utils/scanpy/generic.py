from anndata import AnnData
import scanpy as sc


def _generic_preprocess(
    ad: AnnData, min_cells: int = None, min_genes: int = None, max_genes: int = None, max_pct_mito: int = None
):
    ad.var_names_make_unique()
    sc.pp.calculate_qc_metrics(ad, inplace=True)
    if min_genes:
        sc.pp.filter_cells(ad, min_genes=min_genes)
    if min_cells:
        sc.pp.filter_genes(ad, min_cells=min_cells)
    #     ad = ad[ad.obs.pct_counts_mt < 15]
    if max_genes:
        ad = ad[ad.obs.n_genes_by_counts < max_genes]
    sc.pp.normalize_total(ad, target_sum=1e4)
    return ad
