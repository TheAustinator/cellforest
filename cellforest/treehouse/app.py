from typing import TYPE_CHECKING

from dataforest.treehouse.components.scatter import Scatter, DropDownColumn, SliderColumn, DropDownOptionsStore, \
    ColumnAdder, ScatterAnnotator
from jupyter_dash import JupyterDash

from cellforest.treehouse.Model import CellModel

if TYPE_CHECKING:
    from dataforest.core.DataTree import DataTree


STYLESHEET = ["https://codepen.io/chriddyp/pen/bWLwgP.css"]


def run(tree: "DataTree") -> JupyterDash:
    JupyterDash.infer_jupyter_proxy_config()
    app = JupyterDash(__name__, external_stylesheets=STYLESHEET)
    model = CellModel(tree)
    scatter = Scatter(app)
    dim_dropdowns = DropDownColumn(app, data=dict(
        title_list=model.DIMENSIONS,
        options=model.col_options
    ))
    sliders = SliderColumn(data={"slider_params": model.MARKER_ARGS_PARAMS})
    col_adder = ColumnAdder()
    annotator = ScatterAnnotator(data={"options": model.col_options})
    colname_store = DropDownOptionsStore()

    colname_store.callbacks(annotator, subscribers=[*dim_dropdowns.contents, annotator.dropdown])
    annotator.callbacks(scatter, col_adder)
    scatter.callbacks(dim_dropdowns, sliders)

    return app
