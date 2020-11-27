from dash.dependencies import Input, Output, State
import dash_core_components as dcc
import dash_html_components as html

from dataforest.treehouse.components.core import ModelBlock, Table
from dataforest.treehouse.components.scatter import Scatter, ScatterTab


class GeneData(ModelBlock):
    def layout(self):
        return html.Div(
            className="row",
            children=[
                html.Div(id="top-genes", className="three columns", children=[html.Div(id="top-genes-table")]),
                html.Div(
                    id="dif-genes",
                    className="three columns",
                    children=[
                        dcc.Input(id="dif-genes-min", type="number"),
                        dcc.Input(id="dif-genes-table", type="number"),
                    ],
                ),
            ],
        )

    def callbacks(self, scatter: "Scatter"):
        @self.app.callback(
            Output("top-genes-table", "children"), Input("umap", "selectedData"),
        )
        def show_top_genes(selected_data):
            cells = scatter.get_selected_indices(selected_data)
            df = self.model.get_expr_df(self.rna, cells).round(4)
            self.top_genes = df
            return [Table(data={"df": df})]

        @self.app.callback(
            Output("dif-genes-table", "children"), Input("umap", "selectedData"),
        )
        def show_crude_markers(selected_data):
            cells_sel = self.model.get_cell_ids(selected_data)
            cells_all = self.meta.index
            df = self.model.get_crude_markers(cells_sel, cells_all)
            return [Table(data={"df": df.round(4)})]


class CellScatterTab(ScatterTab):
    COMPONENTS = ScatterTab.COMPONENTS + [GeneData]
