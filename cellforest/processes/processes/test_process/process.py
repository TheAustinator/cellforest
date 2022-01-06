from dataforest.hooks import dataprocess


@dataprocess(plots=True)
def test_process(branch: "CellBranch", run_name: str):
    output_pickle_path = branch[run_name].path_map["rna"]
    open(str(output_pickle_path), "a").close()
