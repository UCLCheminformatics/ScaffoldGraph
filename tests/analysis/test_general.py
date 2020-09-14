"""
scaffoldgraph tests.analysis.test_general
"""

from scaffoldgraph.analysis import get_singleton_scaffolds, get_virtual_scaffolds
from ..test_network import long_test_network


def test_get_virtual_scaffolds(network):
    v = get_virtual_scaffolds(network)
    assert len(v) == 19


def test_get_singleton_scaffolds(network):
    s = get_singleton_scaffolds(network)
    assert len(s) == 3
