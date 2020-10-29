from pathlib import Path

from dataforest.utils.loaders.config import get_config_loader, get_config_options
from dataforest.utils.loaders.update_config import get_config_updater

from cellforest import config as _config_module

_CONFIG_DIR = Path(_config_module.__file__).parent
CONFIG_OPTIONS = get_config_options(_CONFIG_DIR)
load_config = get_config_loader(CONFIG_OPTIONS)

update_config = get_config_updater(load_config)
