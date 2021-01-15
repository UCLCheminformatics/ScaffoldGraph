"""
scaffoldgraph.analysis

The analysis package contains functions for analyzing ScaffoldGraphs
"""

from . import diversity
from . import enrichment
from . import general
from . import representation

from .general import *

__all__ = [
    'get_singleton_scaffolds',
    'get_virtual_scaffolds'
]
