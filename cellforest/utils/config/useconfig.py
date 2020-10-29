from functools import wraps
from typing import AnyStr, Dict, Union, Callable, Any

import pytest
from dataforest import get_current_config
from decorator import decorator

from cellforest.utils.config.functions import update_config


# noinspection PyPep8Naming
class useconfig:
    # TODO: can move to dataforest if can figure out how to move `update_config`
    # TODO: currently won't work in cases where e.g. branch["root"].plot_map is already loaded
    #   -- need to get branch to do branch.copy()
    def __init__(self, config: Union[Dict, AnyStr], persist: bool = False):
        """
        Args:
            config: config to use (arg to cellforest.update_config)
            persist: whether or not to keep updated config after execution
        """
        self._config = config
        self._persist = persist

    def __call__(self, func: Callable) -> Callable:
        def wrapper(_, *args, **kwargs) -> Any:
            prev_config = get_current_config()
            update_config(self._config)
            output = func(*args, **kwargs)
            if not self._persist:
                update_config(prev_config)
            return output

        return decorator(wrapper, func)
