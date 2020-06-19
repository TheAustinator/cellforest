from dataforest.templates.PlotMethods import PlotMethods
import matplotlib.pyplot as plt


class PlotMethodsSC(PlotMethods):
    @staticmethod
    def umap(forest, ax=None, labels="cluster_id", save=False, **kwargs):
        ax = ax if ax else plt.gca()
        meta = forest.meta
        for (name, df) in meta.groupby(labels):
            ax.scatter(df["UMAP_1"], df["UMAP_2"], label=name, s=0.1, alpha=0.1)
        ax.legend()
