"""
scaffoldgraph.analysis
"""

from .representation import calc_average_pairwise_similarity, get_over_represented_scaffold_classes
from .general import get_virtual_scaffolds, get_singleton_scaffolds
from .cse import CSE

__all__ = [
    'calc_average_pairwise_similarity',
    'get_over_represented_scaffold_classes',
    'get_virtual_scaffolds',
    'get_singleton_scaffolds',
    'CSE',
]
