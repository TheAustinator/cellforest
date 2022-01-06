from pathlib import Path

from cellforest.hooks import hooks

__all__ = []

"""
systematically import all hooks from submodules
"""
hooks_files = list(filter(lambda p: p.name == "hooks.py", Path(hooks.__path__[0]).rglob("*.py")))
for path in hooks_files:
    path = "cellforest" + str(path).split("cellforest/cellforest")[1][:-3].replace("/", ".")
    module = __import__(path, fromlist=["*"])
    namespace = [x for x in vars(module) if x.startswith("hook_")]
    if "__all__" in module.__dict__:
        names = module.__dict__["__all__"]
    else:
        # otherwise we import all names that don't begin with _
        names = [x for x in module.__dict__ if not x.startswith("_")]

    # now drag them in
    globals().update({k: getattr(module, k) for k in names})
    __all__ += namespace
pass
