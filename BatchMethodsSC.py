from typing import List, Union

import pandas as pd

from dataforest.BatchMethods import BatchMethods
from gsea.GSEAGroup import GSEAGroup
from cellforest.CellORM import CellForest


class BatchMethodsSC(BatchMethods):
    @staticmethod
    def gsea_bulk(
        orm: CellForest, batch_vars: Union[str, list, set, tuple], overwrite: bool = False, **kwargs
    ) -> pd.DataFrame:
        """
        Run multiple GSEAs over a set of varying conditions as specified by
        `batch_vars`.
        Example:
            GSEA between healthy and diseased over `batch_vars`:
            `{experiment_name, cluster_id}`. In this case, `orm`
            would not be used to subset by either of the `batch_vars`, as this
            is done automatically, but the `orm` would need to include
            `disease_state` for `partition`.

        Args:
            orm:
            batch_vars:
            overwrite:
            **kwargs:

        Returns:
            gsea_grp:
        """
        # orm = orm.at("gsea_bulk")
        gsea_grp = GSEAGroup(orm, batch_vars, "gsea_bulk")
        for grp_vals in gsea_grp.groups:
            BatchMethodsSC.logger.info(f"Running GSEA for {grp_vals}")
            gsea_grp.run_grp(grp_vals, **kwargs)
        return gsea_grp.group_results_df

    @staticmethod
    def gsea_bulk_repeat(
        orm: CellForest, batch_vars: Union[str, list, set, tuple], overwrite: bool = True, n_repeat: int = 20, **kwargs
    ) -> List[pd.DataFrame]:
        group_results_df_list = []
        for i in range(n_repeat):
            if kwargs.get("shuffling", False) or i > 0:
                kwargs["no_plot"] = True
            group_results_df_list.append(BatchMethodsSC.gsea_bulk(orm, batch_vars, overwrite, **kwargs))
        return group_results_df_list
