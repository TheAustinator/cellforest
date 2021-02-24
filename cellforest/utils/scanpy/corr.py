from dataforest.utils.analysis.pairwise import pairwise_metric, pairwise_metric_heat
import matplotlib.pyplot as plt
import numpy as np
from seaborn import heatmap, clustermap
from scipy.spatial.distance import cdist
from scipy.stats import pearsonr
from sklearn.preprocessing import scale

from cellforest.utils.scanpy.gene import get_feature
from cellforest.utils.scanpy.group import groupby


def corr_matrix(ad, features):
    pearson_r = lambda x, y: pearsonr(x, y)[0]
    pearson_p = lambda x, y: pearsonr(x, y)[1]

    feature_dict = {f: get_feature(ad, f) for f in features}
    matrix_r = pairwise_metric(feature_dict, pearson_r)
    matrix_p = pairwise_metric(feature_dict, pearson_p)
    matrix_r = matrix_r.loc[features, features]
    matrix_p = matrix_p.loc[features, features]
    return matrix_r, matrix_p


def corr_heatmap(
    ad, query, batch_var=None, ctrl_genes=True, vmin=-0.5, vmax=0.5, pval=0.05, cluster=False, n_test_mult=1
):
    if ctrl_genes:
        if isinstance(ctrl_genes, bool):
            ctrl_genes = get_ctrl_genes(ad, query)
        query = np.unique([*query, *ctrl_genes])
    if batch_var:
        corr_r_list = list()
        corr_p_list = list()
        test_count = len(ad.obs[batch_var].unique())
        for name, _ad in groupby(ad, batch_var):
            test_count_mult = n_test_mult * test_count
            cm, corr_r, corr_p = corr_heatmap(_ad, query, None, ctrl_genes, vmin, vmax, pval, cluster, test_count_mult)
            cm.ax_heatmap.set_title(name)
            corr_r_list.append(corr_r)
            corr_p_list.append(corr_p)
        corr_r = sum(corr_r_list) / test_count
        corr_p = sum(corr_p_list) / test_count
        cm, _, _ = _corr_heatmap(query, corr_r, corr_p, vmin, vmax, pval, cluster, n_test_mult)
        cm.ax_heatmap.set_title("mean")
        return cm, corr_r, corr_p
    corr_r, corr_p = corr_matrix(ad, query)
    corr_p = corr_p.fillna(np.inf)
    cm, _, _ = _corr_heatmap(query, corr_r, corr_p, vmin, vmax, pval, cluster, n_test_mult)
    return cm, corr_r, corr_p


def _corr_heatmap(query, corr_r, corr_p, vmin, vmax, pval, cluster, n_test_mult=1):
    n_tests = n_test_mult * len(query) ** 2 / 2 - len(query)
    corr_mask = corr_p > pval / n_tests
    corr_r = corr_r.fillna(0)
    plot_func = clustermap if cluster else heatmap
    _ = plt.subplots() if not cluster else None
    return plot_func(corr_r, vmin=vmin, vmax=vmax, mask=corr_mask), corr_r, corr_p


def get_ctrl_genes(ad, query: list, plot=True, n=3, figsize=(10, 5), **kwargs):
    var_mean = np.array(ad.X.mean(axis=0))[0]
    var_vars = np.var(ad.X.toarray(), axis=0)
    # std scale the stats themselves so when calculating distance, each stat gets equity
    var_stat = scale(np.array(np.vstack((var_mean, var_vars))).T)
    gene_inds = np.nonzero(ad.var_names.isin(query))[0]
    gene_stat = var_stat[gene_inds]
    neighbors_dict = dict()
    ctrl_inds = []
    for arr in cdist(gene_stat, var_stat):
        inds = arr.argsort()[: n + 1]
        gene = ad.var_names[inds[0]]
        inds = inds[1:]
        ctrl_inds += inds.tolist()
        neighbors_dict[gene] = ad.var_names[inds].tolist()
    if plot:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize, **kwargs)
        ax1.set_title("probe genes")
        ax2.set_title("control genes")
        plot_gene_stats(var_stat, gene_inds, ax1, **kwargs)
        plot_gene_stats(var_stat, ctrl_inds, ax2, **kwargs)
    ctrl_genes = [g for l in neighbors_dict.values() for g in l]
    return ctrl_genes


def plot_gene_stats(var_stat, gene_inds, ax, s=4, alpha=0.5, **kwargs):
    ax.scatter(var_stat[:, 0], var_stat[:, 1], s=s, alpha=alpha, **kwargs)
    ax.scatter(var_stat[gene_inds, 0], var_stat[gene_inds, 1], s=s, alpha=alpha, **kwargs)
    ax.set_xlabel("normalized mean")
    ax.set_ylabel("normalized var")
    return ax
