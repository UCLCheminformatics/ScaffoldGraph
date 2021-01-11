"""
scaffoldgraph tests.utils.test_misc
"""

import scaffoldgraph as sg

from scaffoldgraph.utils import summary
from . import mock_sdf


SUMMARY_GRAPH = """Type: ScaffoldNetwork
Number of molecule nodes: 2
Number of scaffold nodes: 8
Number of edges: 12
Max hierarchy: 3
Min hierarchy: 1
"""


SUMMARY_NODE = """Node c1ccccc1 has the following properties:
Type: scaffold
Hierarchy: 1
Degree: 2
Parent scaffolds: 
Child scaffolds: O=C1CN=C(c2ccccc2)C=CN1 O=C1CN=Cc2ccccc2N1 O=C1CN=C(c2ccccc2)c2ccccc2N1 O=C1CN=C(c2ccccc2)c2ccsc2N1
Child molecules:
"""


def test_bipartite(sdf_file):
    network = sg.ScaffoldNetwork.from_sdf(sdf_file)
    assert summary(network).strip() == SUMMARY_GRAPH.strip()
    assert summary(network, 'c1ccccc1').strip() == SUMMARY_NODE.strip()
