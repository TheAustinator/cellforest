def check_root_files(branch):
    assert (branch["root"].path / "meta.tsv").exists()
    assert (branch["root"].path / "rna.pickle").exists()
    assert (branch["root"].path / "rna.rds").exists()
