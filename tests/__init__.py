"""
scaffoldgraph tests
"""

import os

import pytest
from rdkit import rdBase

rdBase.DisableLog('rdApp.*')


def test_root_dir():
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(name='sdf_file')
def mock_sdf():
    return


@pytest.fixture(name='smiles_file')
def mock_smiles_file():
    return
