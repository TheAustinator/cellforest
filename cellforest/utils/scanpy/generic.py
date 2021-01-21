from anndata import AnnData
import scanpy as sc


def _generic_preprocess(
    ad: AnnData, min_cells: int = 10, min_genes: int = 200, max_genes: int = 2500, max_pct_mito: int = 30
):
    ad.var_names_make_unique()
    sc.pp.calculate_qc_metrics(ad, inplace=True)
    sc.pp.filter_cells(ad, min_genes=min_genes)
    sc.pp.filter_genes(ad, min_cells=min_cells)
    #     ad = ad[ad.obs.pct_counts_mt < 15]
    ad = ad[ad.obs.n_genes_by_counts < max_genes]
    sc.pp.normalize_total(ad, target_sum=1e4)
    return ad
