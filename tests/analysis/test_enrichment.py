"""
scaffoldgraph tests.analysis.test_enrichment
"""

import pytest
import networkx as nx
import random

from scaffoldgraph.analysis import calc_scaffold_enrichment, compound_set_enrichment
from ..test_network import long_test_network


def test_enrichment(network):

    ks_data = {}
    for molecule in network.get_molecule_nodes():
        ks_data[molecule] = {'activity': random.random()}
    nx.set_node_attributes(network, ks_data)
    enrichment = calc_scaffold_enrichment(network, 'activity')
    entry = list(enrichment.items())[0]
    assert entry[0] in network
    assert 'pval' in entry[1]
    assert 'dmax' in entry[1]
    assert '_total' in entry[1]
    assert type(entry[1]['dmax']) == float
    assert type(entry[1]['_total']) == int
    compound_set_enrichment(network, 'activity', mode='ks')

    binom_data = {}
    for molecule in network.get_molecule_nodes():
        binom_data[molecule] = {'activity': random.choice([0, 1])}
    nx.set_node_attributes(network, binom_data)
    enrichment = calc_scaffold_enrichment(network, 'activity', mode='b')
    entry = list(enrichment.items())[0]
    assert entry[0] in network
    assert 'pval' in entry[1]
    assert '_active' in entry[1]
    assert '_total' in entry[1]
    assert type(entry[1]['_active']) == int
    assert type(entry[1]['_total']) == int
    compound_set_enrichment(network, 'activity', mode='b')

    with pytest.raises(ValueError):
        compound_set_enrichment(network, 'activity', mode='not a mode')
