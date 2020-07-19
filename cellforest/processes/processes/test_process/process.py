from dataforest.hooks import dataprocess


@dataprocess()
def test_process(forest: "CellForest", run_name: str):
    output_pickle_path = forest[run_name].path_map["rna"]
    open(str(output_pickle_path), "a").close()
