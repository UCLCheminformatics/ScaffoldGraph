"""
scaffoldgraph tests.core.test_scaffold
"""

import pytest
import pickle

from rdkit import Chem

from scaffoldgraph.core.scaffold import *


@pytest.fixture(name='scaffold')
def basic_scaffold():
    # murcko scaffold smiles
    mol = Chem.MolFromSmiles('O=C(Nc1ccccc1)Nc1cc2c(s1)CNCC2')
    scaffold = Scaffold(mol)
    return scaffold


def test_new():
    scaffold = Scaffold(None)
    assert scaffold is None


def test_pickle(scaffold):
    b = pickle.dumps(scaffold)
    s = pickle.loads(b)
    assert s.atoms
    assert s.bonds
    assert s.rings
    assert s.ring_systems
    assert s.smiles


def test_smiles(scaffold):
    assert scaffold.smiles == 'O=C(Nc1ccccc1)Nc1cc2c(s1)CNCC2'
    assert scaffold.get_canonical_identifier() == scaffold.smiles
    assert scaffold == Scaffold(Chem.MolFromSmiles(scaffold.smiles))
    assert scaffold == scaffold.smiles
    assert str(scaffold) == scaffold.smiles
    assert hash(scaffold) == hash(scaffold.smiles)


def test_name(scaffold):
    assert scaffold.name is None
    scaffold.name = 'TEST'
    assert scaffold.name == 'TEST'
    assert repr(scaffold) == '<Scaffold at {}>'.format(hex(id(scaffold)))
    assert bool(scaffold) is True


def test_atoms(scaffold):
    atoms = scaffold.atoms
    assert len(atoms) == scaffold.mol.GetNumAtoms()
    assert all([isinstance(x, Chem.Atom) for x in atoms])


def test_bonds(scaffold):
    bonds = scaffold.bonds
    assert len(bonds) == scaffold.mol.GetNumBonds()
    assert all([isinstance(x, Chem.Bond) for x in bonds])


def test_rings(scaffold):
    rings = scaffold.rings
    assert isinstance(rings, RingStack)
    assert hasattr(rings, 'owner')
    assert hasattr(rings, 'info')
    assert hasattr(rings, 'atom_rings')
    assert hasattr(rings, 'bond_rings')
    assert rings.count == 3 and len(rings) == 3
    assert repr(rings) == '<RingStack at {}>'.format(hex(id(rings)))
    assert isinstance(rings[0], Ring)
    assert len([x for x in rings]) == 3
    assert isinstance(rings.info, Chem.RingInfo)
    assert len(rings.atom_rings) == 3 and len(rings.bond_rings) == 3
    ring = rings[1]
    assert hasattr(ring, 'owner')
    assert hasattr(ring, 'aix')
    assert hasattr(ring, 'bix')
    assert all([isinstance(x, Chem.Bond) for x in ring.bonds])
    assert all([isinstance(x, Chem.Atom) for x in ring.atoms])
    assert isinstance(ring.size, int)
    assert len(ring) == len(ring.atoms)
    assert repr(ring) == '<Ring at {}>'.format(hex(id(ring)))
    assert len(ring.get_attachment_points()) == 1
    assert ring.is_exocyclic_attachment(ring.atoms[0]) is False
    assert ring.get_ring_system().size == 9
    assert len(rings.to_list()) == 3
    subset = rings[0:2]
    assert len(subset) == 2
    assert subset[0] != subset[1]


def test_ring_systems(scaffold):
    rings = scaffold.ring_systems
    assert isinstance(rings, RingSystemStack)
    assert hasattr(rings, 'owner')
    assert hasattr(rings, 'ring_indexes')
    assert hasattr(rings, 'atom_rings')
    assert hasattr(rings, 'bond_rings')
    assert rings.count == 2 and len(rings) == 2
    assert repr(rings) == '<RingSystemStack at {}>'.format(hex(id(rings)))
    assert isinstance(rings[0], RingSystem)
    assert len([x for x in rings]) == 2
    assert len(rings.atom_rings) == 2 and len(rings.bond_rings) == 2
    ring = rings[1]
    assert hasattr(ring, 'owner')
    assert hasattr(ring, 'aix')
    assert hasattr(ring, 'bix')
    assert hasattr(ring, 'rix')
    assert all([isinstance(x, Chem.Bond) for x in ring.bonds])
    assert all([isinstance(x, Chem.Atom) for x in ring.atoms])
    assert isinstance(ring.size, int)
    assert len(ring) == len(ring.atoms)
    assert repr(ring) == '<RingSystem at {}>'.format(hex(id(ring)))
    assert isinstance(ring[0], Ring)
    assert len(list(ring.get_rings())) == 2
    assert len(ring.get_attachment_points()) == 1
    assert ring.is_exocyclic_attachment(ring.atoms[0]) is False
    subset = rings[1:]
    assert len(subset) == 1
    assert isinstance(subset[0][0], Ring)
    assert len(subset[0][0:2]) == 2

