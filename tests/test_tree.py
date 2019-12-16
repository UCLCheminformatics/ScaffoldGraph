"""
scaffoldgraph tests.test_tree
"""

import networkx as nx
import pytest

import scaffoldgraph as sg

from . import mock_sdf, mock_smiles_file


@pytest.fixture(name='test_tree')
def test_tree_graph(sdf_file):
    tree = sg.ScaffoldTree.from_sdf(sdf_file)
    return tree


def test_tree_from_sdf(sdf_file):
    tree = sg.ScaffoldTree.from_sdf(sdf_file)
    assert tree.num_scaffold_nodes == 5
    assert tree.num_molecule_nodes == 2
    assert nx.is_tree(tree)


def test_tree_from_smiles(smiles_file):
    tree = sg.ScaffoldTree.from_smiles_file(smiles_file)
    assert tree.num_scaffold_nodes == 5
    assert tree.num_molecule_nodes == 2
    assert nx.is_tree(tree)


def test_repr(test_tree):
    assert repr(test_tree) == '<ScaffoldTree at {}>'.format(hex(id(test_tree)))
