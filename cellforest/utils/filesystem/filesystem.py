import gzip
import os
import shutil

from pathlib import Path


def compress(*args, keep_orig=True):
    for arg in args:
        with open(arg, "rb") as f_in:
            with gzip.open(str(arg) + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        if not keep_orig:
            os.remove(arg)


def compress_move(filenames, src_dir, dst_dir, keep_orig=True):
    src_dir = Path(src_dir)
    dst_dir = Path(dst_dir)
    for f in filenames:
        with open(src_dir / f, "rb") as f_in:
            with gzip.open(str(dst_dir / f) + ".gz", "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        if not keep_orig:
            os.remove(src_dir / f)


def attempt_move(source, dst, filenames, keep=True):
    mover = shutil.copy if keep else shutil.move
    try:
        [mover(str(source / f), str(dst) + "/") for f in filenames]
    except shutil.Error:  # output already exists
        pass
