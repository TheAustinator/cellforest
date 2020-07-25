import os
import shutil
import tarfile
import urllib.request

import pandas as pd
from pathlib import Path

from cellforest import Counts
from cellforest.utils import compress_move

DATA_URL = "https://s3-us-west-2.amazonaws.com/10x.files/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz"


def get_test_data():
    """
    Get sample data from 10X for testing. The data comes in the format of v2
    chemistry, and a v3 version is artificially created, as well as v3 .gz
    version
    Returns:

    """
    _get_test_data_slice(300, 150)


def _get_test_data_slice(n_cells, n_genes):
    # create sample metadata
    data_dir = Path(__file__).parent.parent / "data"
    os.makedirs(data_dir, exist_ok=True)
    subdirs = ["v3_gz/sample_1", "v3_gz/sample_2"]
    sample_metadata = pd.DataFrame(
        {"sample": ["sample_1", "sample_2"], "path_rna": [str(data_dir / x) for x in subdirs]}
    )
    sample_metadata.to_csv(data_dir / "sample_metadata.tsv", sep="\t", index=False)

    # pull and unzip data from 10X
    data_dir_gzip = data_dir / "v3_gz"
    download_data(data_dir)

    # cut data into two samples of 200 cells x 100 genes
    src = data_dir / "filtered_gene_bc_matrices/hg19/"
    files = os.listdir(src)
    rna = Counts.from_cellranger(src)
    rna_1 = rna[:n_cells, :n_genes]
    rna_2 = rna[n_cells : 2 * n_cells, :n_genes]

    # save a v2 chemistry version
    dst_1_v2 = data_dir / "v2/sample_1/"
    dst_2_v2 = data_dir / "v2/sample_2/"
    os.makedirs(dst_1_v2, exist_ok=True)
    os.makedirs(dst_2_v2, exist_ok=True)
    rna_1.to_cellranger(dst_1_v2, gz=False, chemistry="v2")
    rna_2.to_cellranger(dst_2_v2, gz=False, chemistry="v2")

    # save a v3 chemistry version (features.tsv with third column)
    files.remove("genes.tsv")
    files.append("features.tsv")
    dst_1_v3 = data_dir / "v3/sample_1/"
    dst_2_v3 = data_dir / "v3/sample_2/"
    os.makedirs(dst_1_v3, exist_ok=True)
    os.makedirs(dst_2_v3, exist_ok=True)
    rna_1.to_cellranger(dst_1_v3, gz=False, chemistry="v3")
    rna_2.to_cellranger(dst_2_v3, gz=False, chemistry="v3")

    # save a gzipped version
    dst_1_gz = data_dir_gzip / "sample_1/"
    dst_2_gz = data_dir_gzip / "sample_2/"
    os.makedirs(dst_1_gz, exist_ok=True)
    os.makedirs(dst_2_gz, exist_ok=True)
    compress_move(files, dst_1_v3, dst_1_gz)
    compress_move(files, dst_2_v3, dst_2_gz)

    # remove downloads
    shutil.rmtree(data_dir / "filtered_gene_bc_matrices",)
    os.remove(data_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz")


def download_data(data_dir):
    download_path = data_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz"
    urllib.request.urlretrieve(
        DATA_URL, filename=download_path,
    )
    tar = tarfile.open(download_path, "r:gz")
    tar.extractall(data_dir)
    tar.close()


def get_test_data_full():
    """
    Similar to `get_test_data`, but gets full dataset without slicing and saves
    it to data/full
    """
    data_dir = Path(__file__).parent.parent / "data"
    download_data(data_dir)
    src = data_dir / "filtered_gene_bc_matrices/hg19/"
    dst = data_dir / "full"
    dst.mkdir(exist_ok=True)
    rna = Counts.from_cellranger(src)
    rna.to_cellranger(dst, gz=False, chemistry="v3")
    shutil.rmtree(data_dir / "filtered_gene_bc_matrices",)
    os.remove(data_dir / "pbmc3k_filtered_gene_bc_matrices.tar.gz")


if __name__ == "__main__":
    get_test_data()
