"""
scaffoldgraph tests.utils.test_aggregate
"""

import scaffoldgraph as sg

from scaffoldgraph.utils import aggregate
from .. import mock_sdf, mock_sdf_2


def test_aggregate(sdf_file, sdf_file_2):
    net_1 = sg.ScaffoldNetwork.from_sdf(sdf_file)
    net_2 = sg.ScaffoldNetwork.from_sdf(sdf_file_2)
    network = aggregate([net_1, net_2])
    assert network.num_scaffold_nodes == 14
    assert network.num_molecule_nodes == 4
