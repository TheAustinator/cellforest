import os

from tests.test_init import *
from tests.test_rna import *
from tests.test_process import *

from path import Path  # pathlib doesn't have rmtree

keep_dirs = ["v2", "v3", "v3_gz"]
data_dirs = filter(lambda x: x.isdir(), list(Path("data").glob("*")))
data_dirs = filter(lambda x: x.name not in keep_dirs, data_dirs)
list(map(lambda x: x.rmtree(), data_dirs))
