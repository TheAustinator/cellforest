import os
import pickle

from pathlib import Path
from scipy.io import mmread

from cellforest.r_utils.shell_command import shell_command


RDS_CONVERTER_SCRIPT = Path(__file__).parent.absolute().parent / "r_scripts" / "rds_converter.R"


def gzip_replace(filepath):
    import gzip
    import shutil
    filepath = str(filepath)
    with open(filepath, 'rb') as f_in:
        with gzip.open(filepath + '.gz', 'wb') as f_out:
            shutil.copyfileobj(f_in, f_out)
    os.remove(filepath)


def seurat_rds_to_sparse_pickle(rds_path, output_dir):
    output_dir = Path(output_dir)
    if not os.path.isfile(rds_path):
        raise ValueError(f"rds file does not exists: {rds_path}")
    output_mtx_path = output_dir / "matrix.mtx"
    output_pickle_path = output_dir / "matrix.pickle"
    output_cell_ids_path = output_dir / "cell_ids.tsv"
    output_genes_path = output_dir / "genes.tsv"
    output_metadata_path = output_dir / "cell_metadata.tsv"
    cmd_str = f"Rscript {RDS_CONVERTER_SCRIPT} {rds_path} {output_mtx_path} {output_cell_ids_path} {output_genes_path} {output_metadata_path}"
    shell_command(cmd_str, output_dir, "matrix_rds_to_sparse_pickle")
    sparse_matrix = mmread(str(output_mtx_path)).T.tocsr()
    with open(str(output_pickle_path), "wb") as f:
        pickle.dump(sparse_matrix, f, protocol=pickle.HIGHEST_PROTOCOL)
    paths_to_zip = [output_mtx_path, output_cell_ids_path, output_genes_path]
    for filepath in paths_to_zip:
        gzip_replace(filepath)
