# your_package/__init__.py
# scCNA/__init__.py

from .CNAfinder import find_cnas
from .utils import annotate_gene_position
from .goldstandard import simulate_cnas_windowed

__version__ = "0.1.0"
