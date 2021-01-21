import matplotlib.pyplot as plt

from cellforest import CellBranch


def umi_hist_facet_stratify(br: CellBranch, facet: str, stratify: str, histtype="step", xlim=(0, 12500), **kwargs):
    # TODO: use subplots outside of loop and return single ax
    for indication, brr in br.groupby(facet):
        print(indication)
        fig, ax = plt.subplots(figsize=(20, 5))
        for lane, brrr in brr.groupby(stratify):
            brrr.rna.hist(bins=200, ax=ax, labels=brrr.meta[stratify], histtype=histtype, **kwargs)
        ax.set_xlim(*xlim)
        plt.show()
