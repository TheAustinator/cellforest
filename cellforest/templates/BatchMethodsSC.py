from typing import List, Union

from dataforest.processes.core.BatchMethods import BatchMethods

# from scgsea.GSEAGroup import GSEAGroup
import pandas as pd

from cellforest.templates.CellBranch import CellBranch


class BatchMethodsSC(BatchMethods):
    @staticmethod
    def gsea_bulk(
        forest: CellBranch, batch_vars: Union[str, list, set, tuple], overwrite: bool = False, **kwargs,
    ) -> pd.DataFrame:
        """
        Run multiple GSEAs over a set of varying conditions as specified by
        `batch_vars`.
        Example:
            GSEA between healthy and diseased over `batch_vars`:
            `{experiment_name, cluster_id}`. In this case, `branch`
            would not be used to subset by either of the `batch_vars`, as this
            is done automatically, but the `branch` would need to include
            `disease_state` for `partition`.

        Args:
            forest:
            batch_vars:
            overwrite:
            **kwargs:

        Returns:
            gsea_grp:
        """
        gsea_grp = GSEAGroup(forest, batch_vars, "gsea")
        for grp_vals in gsea_grp.groups:
            BatchMethodsSC.logger.info(f"Running GSEA for {grp_vals}")
            gsea_grp.run_grp(grp_vals, **kwargs)
        return gsea_grp.group_results_df

    @staticmethod
    def gsea_bulk_repeat(
        forest: CellBranch,
        batch_vars: Union[str, list, set, tuple],
        overwrite: bool = True,
        n_repeat: int = 20,
        **kwargs,
    ) -> List[pd.DataFrame]:
        group_results_df_list = []
        for i in range(n_repeat):
            if kwargs.get("shuffling", False) or i > 0:
                kwargs["no_plot"] = True
            group_results_df_list.append(BatchMethodsSC.gsea_bulk(forest, batch_vars, overwrite, **kwargs))
        return group_results_df_list
