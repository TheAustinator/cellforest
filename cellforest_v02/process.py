def _load_obs(*paths):
    # sequentially load and merge, concat on axis=1, and overwrite
    pass


def _load_mapping():
    pass


def _load_obsm():
    # use _load_mapping for this and other key val stores
    pass


def _load_obsp():
    pass


def _load_X():
    # use for raw and X
    # handle sparse and dense
    pass


def _load_uns():
    pass


class process:
    # even worth having the dataforest abstracted version, or will that just slow down?
    ATTR_LOADERS = {
        "obs": _load_obs,
        "obsm": _load_obsm,
    }
