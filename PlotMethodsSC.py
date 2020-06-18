import matplotlib.pyplot as plt

from dataforest.PlotMethods import PlotMethods


class PlotMethodsSC(PlotMethods):
    @staticmethod
    def umap(orm, ax=None, labels="cluster_id", save=False, **kwargs):
        ax = ax if ax else plt.gca()
        meta = orm.meta
        for (name, df) in meta.groupby(labels):
            ax.scatter(df["UMAP_1"], df["UMAP_2"], label=name, s=0.1, alpha=0.1)
        ax.legend()
