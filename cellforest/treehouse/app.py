from typing import TYPE_CHECKING

from dataforest.treehouse.components.scatter import Scatter, DropDownColumn, SliderColumn, DropDownOptionsStore, \
    ColumnAdder, ScatterAnnotator
from jupyter_dash import JupyterDash
import pandas as pd

from cellforest.treehouse.Model import CellModel
from cellforest.treehouse.components import GeneData

if TYPE_CHECKING:
    from dataforest.core.DataTree import DataTree


STYLESHEET = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]


class TreeHouse:
    def __init__(self, tree: "DataTree"):
        self.tree = tree
        JupyterDash.infer_jupyter_proxy_config()
        self.app = JupyterDash(__name__, external_stylesheets=STYLESHEET)
        self.server = self.app.server
        self.model = CellModel(tree)
        self.scatter = Scatter(self.app, data={"model": self.model})
        self.dim_dropdowns = DropDownColumn(self.app, data=dict(
            title_list=self.model.DIMENSIONS,
            options=self.model.col_options
        ))
        self.sliders = SliderColumn(data={"slider_params": self.model.MARKER_ARGS_PARAMS})
        self.col_adder = ColumnAdder()
        self.annotator = ScatterAnnotator(data={"model": self.model, "options": self.model.col_options})
        self.colname_store = DropDownOptionsStore()
        self.gene_data = GeneData(self.app, data={"model": self.model})

        self.colname_store.callbacks(master=self.annotator, subscribers=[*self.dim_dropdowns.contents, self.annotator.dropdown])
        self.annotator.callbacks(self.scatter, self.col_adder)
        self.scatter.callbacks(self.dim_dropdowns, self.sliders)
        self.gene_data.callbacks(self.scatter)

    @property
    def top_genes(self):
        return self.model.top_genes

    @property
    def diff_genes(self):
        return self.model.diff_genes

    @property
    def selected_cells(self):
        return pd.Series(getattr(self.scatter, "selected_indices", None))
