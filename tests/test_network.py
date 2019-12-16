"""
scaffoldgraph tests.test_network
"""

import pytest

import scaffoldgraph as sg

from . import mock_sdf, mock_smiles_file


@pytest.fixture(name='test_net')
def test_network(sdf_file):
    network = sg.ScaffoldNetwork.from_sdf(sdf_file)
    return network


def test_network_from_sdf(sdf_file):
    network = sg.ScaffoldNetwork.from_sdf(sdf_file)
    assert network.num_scaffold_nodes == 8
    assert network.num_molecule_nodes == 2


def test_network_from_smiles(smiles_file):
    network = sg.ScaffoldNetwork.from_smiles_file(smiles_file)
    assert network.num_scaffold_nodes == 8
    assert network.num_molecule_nodes == 2


def test_hiers(sdf_file):
    network = sg.HierS.from_sdf(sdf_file)
    assert network.num_scaffold_nodes == 5
    assert network.num_molecule_nodes == 2


def test_repr(test_net):
    assert repr(test_net) == '<ScaffoldNetwork at {}>'.format(hex(id(test_net)))
