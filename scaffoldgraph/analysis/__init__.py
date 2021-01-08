"""
scaffoldgraph.analysis

The analysis package contains functions for analyzing ScaffoldGraphs
"""

from .representation import calc_average_pairwise_similarity, get_over_represented_scaffold_classes
from .enrichment import calc_scaffold_enrichment, compound_set_enrichment
from .general import get_virtual_scaffolds, get_singleton_scaffolds
from .diversity import diversity_pick_for_scaffold_class


__all__ = [
    'calc_average_pairwise_similarity',
    'get_over_represented_scaffold_classes',
    'calc_scaffold_enrichment',
    'compound_set_enrichment',
    'get_virtual_scaffolds',
    'get_singleton_scaffolds',
    'diversity_pick_for_scaffold_class',
]
