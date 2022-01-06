from scanpy import read_10x_h5, read_10x_mtx


def read_10x(path):
    try:
        ad = read_10x_h5(path)
    except OSError:
        try:
            ad = read_10x_mtx(str(path)[:-3])
        except OSError:
            raise OSError(f"Neither .h5 or .mtx found for {path}")
    return ad
