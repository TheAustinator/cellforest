from tests.test_init import *
from tests.test_rna import *
from tests.test_process import *
from tests.test_data_ops import *
from tests.test_historical import *
from tests.test_version_control import *
from tests.test_datatree import *
from tests.test_cellforestR import *
from tests.test_qc_plots import *

from path import Path  # pathlib doesn't have rmtree

keep_dirs = ["v2", "v3", "v3_gz"]
data_dirs = filter(lambda x: x.isdir(), list(Path("data").glob("*")))
data_dirs = filter(lambda x: x.name not in keep_dirs, data_dirs)
list(map(lambda x: x.rmtree(), data_dirs))
