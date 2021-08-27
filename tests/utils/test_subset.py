"""
scaffoldgraph tests.utils.test_subset

"""
import random
import pytest

from scaffoldgraph.utils.subset import split_graph_by_molecule_attribute
from collections import defaultdict

from ..test_network import long_test_network


def test_split_by_attribute(network):
    key, attrs = 'ATTR', ['ATTR_1', 'ATTR_2', 'ATTR_3']
    with pytest.raises(ValueError):
        split_graph_by_molecule_attribute(network, True, None)
        split_graph_by_molecule_attribute(network, False, None)
    assigned = defaultdict(int)
    for _, mol_data in network.get_molecule_nodes(True):
        attr = random.choice(attrs)
        mol_data[key] = attr
        assigned[attr] += 1
    subgraphs = split_graph_by_molecule_attribute(network, key, None)
    assert len(subgraphs) == len(assigned)
    for u_attr in assigned.keys():
        subgraph = subgraphs[u_attr]
        assert subgraph.num_molecule_nodes == assigned[u_attr]
        assert subgraph.num_scaffold_nodes > 0
        assert all([d == u_attr for n, d in subgraph.get_molecule_nodes(key)])
