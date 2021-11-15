import os
import sys

import yaml

import initialize as init

logger = init.get_logger(__file__)

_RESLIB_FILEPATH = os.path.join(os.path.dirname(__file__), 'lib', 'reslib.yml')
with open(_RESLIB_FILEPATH) as f:
	_RESLIB = yaml.load(f, Loader=yaml.FullLoader)