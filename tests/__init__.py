"""
scaffoldgraph tests
"""

import os

import pytest
from rdkit import Chem
from rdkit import rdBase

rdBase.DisableLog('rdApp.*')


def test_root_dir():
    return os.path.dirname(os.path.abspath(__file__))


@pytest.fixture(name='sdf_file')
def mock_sdf(tmp_path):
    d = tmp_path / "test_data"
    d.mkdir()
    p = d / "test.sdf"
    writer = Chem.SDWriter(str(p))
    writer.write(Chem.MolFromSmiles('CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3'))
    writer.write(Chem.MolFromSmiles('CCC1=CC2=C(S1)N(C(=O)CN=C2C3=CC=CC=C3Cl)C'))
    writer.close()
    return str(p)


@pytest.fixture(name='sdf_file_2')
def mock_sdf_2(tmp_path):
    d = tmp_path / "test_data"
    try:
        d.mkdir()
    except FileExistsError:
        pass
    p = d / "test_2.sdf"
    writer = Chem.SDWriter(str(p))
    writer.write(Chem.MolFromSmiles('C1C(=O)NC2=C(C=C(C=C2)Br)C(=N1)C3=CC=CC=N3'))
    writer.write(Chem.MolFromSmiles('CC1=NN(C2=C1C(=NCC(=O)N2C)C3=CC=CC=C3F)C'))
    writer.close()
    return str(p)


@pytest.fixture(name='smiles_file')
def mock_smiles_file(tmp_path):
    d = tmp_path / "test_data"
    d.mkdir()
    p = d / "test.smi"
    writer = Chem.SmilesWriter(str(p))
    writer.write(Chem.MolFromSmiles('CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3'))
    writer.write(Chem.MolFromSmiles('CCC1=CC2=C(S1)N(C(=O)CN=C2C3=CC=CC=C3Cl)C'))
    writer.close()
    return str(p)


def canon(smiles):
    """Canonicalize SMILES for safety. If canonicalization ever changes this should remain consistent"""
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
