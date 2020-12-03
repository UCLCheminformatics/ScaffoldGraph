"""
scaffoldgraph tests.core.test_fragment
"""

import pytest
from rdkit import Chem

from scaffoldgraph.core.fragment import *


@pytest.fixture(name='mol')
def test_molecule():
    smiles = 'CCN1CCc2c(C1)sc(NC(=O)Nc3ccc(Cl)cc3)c2C#N'
    return Chem.MolFromSmiles(smiles)


def canon(smiles):
    """Canonicalize SMILES for safety. If canonicalization ever changes this should remain consistent"""
    return Chem.MolToSmiles(Chem.MolFromSmiles(smiles))


def test_murcko(mol):
    murcko = get_murcko_scaffold(mol, generic=False)
    assert Chem.MolToSmiles(murcko) == canon('O=C(Nc1ccccc1)Nc1cc2c(s1)CNCC2')
    murcko = get_murcko_scaffold(mol, generic=True)
    assert Chem.MolToSmiles(murcko) == canon('CC(CC1CCCCC1)CC1CC2CCCCC2C1')


def test_annotation(mol):
    annotation = Chem.MolToSmiles(get_annotated_murcko_scaffold(mol))
    annotation = annotation.replace('1*', '*')
    annotation = annotation.replace('2*', '*')
    annotation = annotation.replace('3*', '*')
    assert annotation.count('*') == 3


def test_murcko_all(mol):
    frags = get_all_murcko_fragments(mol, break_fused_rings=True)
    assert len(frags) == 6
    frags = get_all_murcko_fragments(mol, break_fused_rings=False)
    assert len(frags) == 3


def test_murcko_next(mol):
    scf = get_murcko_scaffold(mol)
    frags_1 = get_next_murcko_fragments(scf, break_fused_rings=True)
    frags_1 = {Chem.MolToSmiles(x) for x in frags_1}
    assert len(frags_1) == 2
    frags_2 = get_next_murcko_fragments(scf, break_fused_rings=False)
    frags_2 = {Chem.MolToSmiles(x) for x in frags_2}
    assert len(frags_2) == 2
    assert len(frags_1.intersection(frags_2)) == 1


def test_collect_linker_atoms():
    mol = Chem.MolFromSmiles('CCCCCCCCCc1ccccc1')
    remove_atoms = set()
    a = collect_linker_atoms(mol.GetAtomWithIdx(0), remove_atoms, True)
    assert len(a) == 1
    assert len(remove_atoms) == 9
    remove_atoms.clear()
    a = collect_linker_atoms(mol.GetAtomWithIdx(0), remove_atoms, False)
    assert len(a) == 1
    assert len(remove_atoms) == 8


def test_remove_exocylic_attachments(mol):
    edited = remove_exocyclic_attachments(mol)
    assert Chem.MolToSmiles(edited) == canon('CCN1CCc2c(sc(NCNc3ccc(Cl)cc3)c2C#N)C1')


def test_genericise_scaffold(mol):
    generic = genericise_scaffold(mol)
    assert Chem.MolToSmiles(generic) == canon('CCC1CCC2C(C1)CC(CC(C)CC1CCC(C)CC1)C2CC')
