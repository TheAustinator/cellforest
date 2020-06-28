import gzip
import os
import shutil
import tarfile
import urllib.request

import pandas as pd
from pathlib import Path


def get_test_data():
    data_dir = Path(__file__).parent.parent / "data"
    data_dir_gzip = data_dir / "v3_gz"
    download_path = data_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz"
    urllib.request.urlretrieve(
        "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz",
        filename=download_path,
    )
    tar = tarfile.open(download_path, "r:gz")
    tar.extractall(data_dir)
    tar.close()
    source = data_dir / "filtered_gene_bc_matrices" / "hg19/"
    dst_1_v2 = data_dir / "v2/sample_1/"
    dst_2_v2 = data_dir / "v2/sample_2/"
    os.makedirs(dst_1_v2, exist_ok=True)
    os.makedirs(dst_2_v2, exist_ok=True)
    files = os.listdir(source)
    _attempt_move(source, dst_1_v2, files)
    _attempt_move(source, dst_2_v2, files)

    # create v3 version of files (features.tsv with third column)
    df = pd.read_csv(dst_1_v2 / "genes.tsv", sep="\t", header=None)
    df["mode"] = "Gene Expression"
    dst_1 = data_dir / "v3/sample_1/"
    dst_2 = data_dir / "v3/sample_2/"
    os.makedirs(dst_1, exist_ok=True)
    os.makedirs(dst_2, exist_ok=True)
    df.to_csv(dst_1 / "features.tsv", header=None, index=False, sep="\t")
    df.to_csv(dst_2 / "features.tsv", header=None, index=False, sep="\t")
    files.remove("genes.tsv")
    _attempt_move(source, dst_1, files)
    _attempt_move(source, dst_2, files)
    files.append("features.tsv")

    # add gzipped versions
    dst_1_gz = data_dir_gzip / "sample_1/"
    dst_2_gz = data_dir_gzip / "sample_2/"
    os.makedirs(dst_1_gz, exist_ok=True)
    os.makedirs(dst_2_gz, exist_ok=True)
    for dst in [dst_1_gz, dst_2_gz]:
        for f in files:
            with open(dst_1 / f, "rb") as f_in:
                with gzip.open(str(dst / f) + ".gz", "wb") as f_out:
                    shutil.copyfileobj(f_in, f_out)

    # remove downloads
    shutil.rmtree(data_dir / "filtered_gene_bc_matrices",)
    os.remove(data_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz")


def _attempt_move(source, dst, filenames, keep=True):
    mover = shutil.copy if keep else shutil.move
    try:
        [mover(str(source / f), str(dst) + "/") for f in filenames]
    except shutil.Error:  # output already exists
        pass


if __name__ == "__main__":
    get_test_data()
