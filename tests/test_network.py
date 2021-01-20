"""
scaffoldgraph tests.test_network
"""

import pytest
import os

from pathlib import Path

import scaffoldgraph as sg

from . import mock_sdf, mock_smiles_file


TEST_DATA_DIR = Path(__file__).resolve().parent / 'data'


@pytest.fixture(name='test_net')
def test_network(sdf_file):
    network = sg.ScaffoldNetwork.from_sdf(sdf_file)
    return network


@pytest.fixture(name='network')
def long_test_network():
    network = sg.ScaffoldNetwork.from_smiles_file(str(TEST_DATA_DIR / 'test_smiles.smi'))
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


def test_hierarchy_functions(network):
    hierarchy_sizes = network.get_hierarchy_sizes()
    assert hierarchy_sizes[1] == 7
    assert hierarchy_sizes[2] == 10
    assert hierarchy_sizes[3] == 7
    assert hierarchy_sizes[4] == 1
    assert network.max_hierarchy() == 4
    assert network.min_hierarchy() == 1
    s_in_h2 = {
        'C1=Cn2cnnc2CN=C1', 'O=C1CN=C(c2ccccn2)C=CN1', 'O=C1CN=Cc2ccccc2N1',
        'C1=CC(c2ccccc2)=[NH+]CC=N1', 'C1=Nc2ccccc2C=[NH+]C1',
        'O=C1CC(=O)N(c2ccccc2)C=CN1', 'O=C1CC(=O)Nc2ccccc2N1',
        'O=C1CN=C(c2ccccc2)C=CN1', 'O=C1C[NH+]=Cc2ccccc2N1',
        'O=C1C[NH+]=C(c2ccccc2)C=CN1'
    }
    assert s_in_h2 == set(network.get_scaffolds_in_hierarchy(2))


def test_simple_functions(network):
    assert network.scaffold_in_graph('C1=Cn2cnnc2CN=C1') is True
    # Below is the non-canonical SMILES of the above
    assert network.scaffold_in_graph('C1=C-n2:c:n:n:c:2-C-N=C-1') is True
    assert network.scaffold_in_graph('c1ccccc1CCNc2ccccc2') is False
    assert network.molecule_in_graph('Adinazolam') is True
    assert network.molecule_in_graph('Citalopram') is False


def test_traversal(network):
    s_for_adinazolam = {
        'c1ccc(C2=NCc3nncn3-c3ccccc32)cc1', 'C1=NCc2nncn2-c2ccccc21',
        'C1=Cn2cnnc2CN=C1c1ccccc1', 'C1=Cn2cnnc2CN=C1', 'c1nnc[nH]1'
    }
    assert set(network.get_scaffolds_for_molecule('Adinazolam')) == s_for_adinazolam
    m_for_scaffold = {'Adinazolam', 'Alprazolam'}
    assert set(network.get_molecules_for_scaffold('c1nnc[nH]1')) == m_for_scaffold


def test_separate_disconnected(network):
    assert len(network.separate_disconnected_components(sort=True)) == 2
    assert type(network.separate_disconnected_components()[0]) == type(network)


def test_repr(test_net):
    assert repr(test_net) == '<ScaffoldNetwork at {}>'.format(hex(id(test_net)))
