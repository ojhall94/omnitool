import os
import inspect
from pkg_resources import resource_filename

__ROOT__ = '/'.join(os.path.abspath(inspect.getfile(inspect.currentframe())).split('/')[:-1])
# default data directory
datadir = resource_filename('omnitool', 'data')
