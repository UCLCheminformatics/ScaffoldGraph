"""
scaffoldgraph tests.utils.test_bipartite
"""

import scaffoldgraph as sg
import networkx as nx

from scaffoldgraph.utils.bipartite import make_bipartite_graph
from . import mock_sdf


def test_bipartite(sdf_file):
    network = sg.ScaffoldNetwork.from_sdf(sdf_file)
    biparite = make_bipartite_graph(network)
    assert nx.is_bipartite(biparite)
