"""
scaffoldgraph
"""

from loguru import logger

from . import prioritization
from . import utils
from .core import get_next_murcko_fragments, get_all_murcko_fragments, get_murcko_scaffold
from .network import ScaffoldNetwork, HierS
from .tree import ScaffoldTree, tree_frags_from_mol

__version__ = '1.0.1'

__all__ = [
    '__version__',
    'HierS',
    'ScaffoldNetwork',
    'ScaffoldTree',
    'tree_frags_from_mol',
    'get_next_murcko_fragments',
    'get_all_murcko_fragments',
    'get_murcko_scaffold',
]

logger.disable(__name__)
