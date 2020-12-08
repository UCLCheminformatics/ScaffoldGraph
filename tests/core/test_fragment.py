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
    murcko = get_murcko_scaffold(mol, generic=True, remove_exocyclic=True)
    assert Chem.MolToSmiles(murcko) == canon('C1CCC(CCCC2CC3CCCCC3C2)CC1')
    murcko = get_murcko_scaffold(mol, generic=True, remove_exocyclic=True, collapse_linkers=True)
    assert Chem.MolToSmiles(murcko) == canon('C1CCC(C2CC3CCCCC3C2)CC1')


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


def test_linker_collapse(mol):
    from scaffoldgraph.core.fragment import _collapse_linker_bonds
    collapsed = _collapse_linker_bonds(mol, retain_het=False)
    assert Chem.MolToSmiles(collapsed) == canon('CN1CCc2c(sc(C(=O)c3ccc(Cl)cc3)c2N)C1')
    collapsed = _collapse_linker_bonds(mol, retain_het=True)
    assert Chem.MolToSmiles(collapsed) == canon('CN1CCc2c(sc(NC(=O)Nc3ccc(Cl)cc3)c2N)C1')


def test_ring_toplogy():
    # Replicate figure 1 from paper: Scaffold Topologies II: Analysis of Chemical Databases
    smiles = 'CC(C)c1ccc(C)cc1OC(=O)C2(CCC3C2)C(=C)C3(C)C'
    mol = Chem.MolFromSmiles(smiles)
    topology = get_ring_toplogy_scaffold(mol)
    assert Chem.MolToSmiles(topology) == canon('C1CC1C12CC1C2')


def _test_topology_helper(smiles, expected):
    mol = Chem.MolFromSmiles(smiles)
    topology = get_ring_toplogy_scaffold(mol)
    assert Chem.MolToSmiles(topology) == canon(expected)


def test_ring_topology_extended():
    # Replicate figure 2 from paper: Scaffold Topologies II: Analysis of Chemical Databases
    # Figure 2a: topologies, Figure 2b: Examples of molecules with each topology
    # First 10 examples
    _test_topology_helper('NCC1(CC(=O)O)CCCCC1', 'C1CC1')  # 1
    _test_topology_helper('CNS(=O)(=O)Cc1ccc2[nH]cc(CCN(C)C)c2c1', 'C1C2CC12')  # 2
    _test_topology_helper('COc1ccc(C(CN(C)C)C2(O)CCCCC2)cc1', 'C1CC1C1CC1')  # 3
    _test_topology_helper('[NH3+][Pt]1([NH3+])OC(=O)C2(CC2)C(=O)O1', 'C1CC12CC2')  # 4
    _test_topology_helper('CC1CCc2cc(F)cc3c(=O)c(C(=O)O)cn1c23', 'C12C3C1C23')  # 5
    _test_topology_helper('NC(=O)N1C2C=CC=CC2CC(=O)C2C=CC=CC21', 'C1C2C1C1CC21')  # 6
    _test_topology_helper('Cc1ccc(-c2cc(C(F)(F)F)nn2-c2ccc(S(C)(=O)=O)cc2)cc1', 'C1CC1C1CC1C1CC1')  # 7
    _test_topology_helper('COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1', 'C1CC1C1C2CC21')  # 8
    _test_topology_helper('O=C(O)COCCN1CCN(C(c2ccccc2)c2ccc(Cl)cc2)CC1', 'C1CC1C(C1CC1)C1CC1')  # 9
    _test_topology_helper('O=C1O[Pt]2(NC3CCCCC3N2)OC1=O', 'C1C2C1C21CC1')  # 10


def _test_connectivity_helper(smiles, expected, single=False):
    mol = Chem.MolFromSmiles(smiles)
    connectivity = get_ring_connectivity_scaffold(mol, single)
    assert Chem.MolToSmiles(connectivity) == canon(expected)


def test_ring_connectivity():
    # Test cases from Figure 1 from the paper: Scaffold analysis of pubchem database
    # as a background for hierarchial scaffold-based visualisation.
    _test_connectivity_helper('CC(C)CC1=CC=C(C=C1)C(C)C(O)=O', '*', False)
    _test_connectivity_helper('CC(C)CC1=CC=C(C=C1)C(C)C(O)=O', '*', True)
    _test_connectivity_helper('CC1=CC(NS(=O)(=O)C2=CC=C(N)C=C2)=NO1', '**', False)
    _test_connectivity_helper('CC1=CC(NS(=O)(=O)C2=CC=C(N)C=C2)=NO1', '**', True)
    _test_connectivity_helper('CN1C2=C(C=C(Cl)C=C2)C(=NCC1=O)C1=CC=CC=C1', '**=*', False)
    _test_connectivity_helper('CN1C2=C(C=C(Cl)C=C2)C(=NCC1=O)C1=CC=CC=C1', '***', True)
    db00741 = '[H][C@@]12CC[C@](O)(C(=O)CO)[C@@]1(C)C[C@H](O)[C@@]1([H])[C@@]2([H])CCC2=CC(=O)CC[c@]12C'
    _test_connectivity_helper(db00741, '*=*=*=*', False)
    _test_connectivity_helper(db00741, '****', True)
