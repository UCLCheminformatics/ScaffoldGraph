"""
scaffoldgraph.core.scaffold
A module defining scaffold objects used in ScaffoldGraph
"""

import weakref

from rdkit.Chem import MolToSmiles


class Scaffold(object):
    """A convenience class for defining scaffolds in scaffoldgraph.
    The class contains methods and properties useful during fragmentation
    and prioritization.

    Scaffold is only used internally by scaffoldgraph
    """

    __slots__ = (
        'mol',
        '_atoms',
        '_bonds',
        '_rings',
        '_ring_systems',
        '_smiles',
        'hash_func',
        '__weakref__',
    )

    def __new__(cls, rdmol, *args, **kwargs):
        if rdmol is None:
            return None
        return super(Scaffold, cls).__new__(cls)

    def __init__(self, rdmol, hash_func=None):
        """Initialize a Scaffold object

        Parameters
        ----------
        rdmol: an rdkit Mol object
        hash_func: a callable function for hash calculation, by default
            the canonical smiles string is used. (optional, default: None)
        """

        self.mol = rdmol
        self._atoms = None
        self._bonds = None
        self._rings = None
        self._ring_systems = None
        self.hash_func = hash_func

    @property
    def name(self):
        """Returns the name/identifier of the scaffold"""
        try:
            return self.mol.GetProp('_Name')
        except KeyError:
            return None

    @name.setter
    def name(self, value):
        """Sets the name/identifier of the scaffold"""
        self.mol.SetProp('_Name', str(value))

    @property
    def atoms(self):
        """Return a list of atoms in the scaffold"""
        if self._atoms is None:
            self._atoms = [a for a in self.mol.GetAtoms()]
        return self._atoms

    @property
    def bonds(self):
        """Return a list of bonds in the scaffold"""
        if self._bonds is None:
            self._bonds = [b for b in self.mol.GetBonds()]
        return self._bonds

    @property
    def rings(self):
        """Return ring information for rings in the scaffold"""
        if self._rings is None:
            self._rings = RingStack(self)
        return self._rings

    @property
    def ring_systems(self):
        """Return ring system information for scaffold"""
        if self._ring_systems is None:
            self._ring_systems = RingSystemStack(self)
        return self._ring_systems

    @property
    def smiles(self):
        """Returns the canonical smiles string of the scaffold"""
        if not hasattr(self, '_smiles'):
            setattr(self, '_smiles', MolToSmiles(self.mol))
        return self._smiles

    @property
    def prioritization_rule(self):
        """Returns the prioritization rule used to select this scaffold"""
        try:
            return self.mol.GetProp('prioritization_rule')
        except KeyError:
            return None

    @prioritization_rule.setter
    def prioritization_rule(self, value):
        """Sets the prioritization rule used to select this scaffold"""
        self.mol.SetProp('prioritization_rule', str(value))

    @property
    def removed_ring_idx(self):
        """Returns the index of the removed ring in the child molecule"""
        try:
            return self.mol.GetIntProp('removed_ring_idx')
        except KeyError:
            return None

    @removed_ring_idx.setter
    def removed_ring_idx(self, value):
        """Sets the index of the removed ring in the child molecule"""
        self.mol.SetIntProp('removed_ring_idx', int(value))

    def get_canonical_identifier(self):
        """Returns a canonical identifier for the scaffold:
        A function for representation calculation can be supplied,
        hash_func = callable
        """
        if self.hash_func:
            return self.hash_func(self.mol)
        return self.smiles

    def __bool__(self):
        """Returns True if the molecule contains at least 1 atom"""
        return len(self.atoms) >= 1

    def __hash__(self):
        """Returns a hash of the canonical identifier. This identifier
        can be supplied in the class initializer"""
        return hash(self.get_canonical_identifier())

    def __eq__(self, other):
        if isinstance(other, str):
            return self.get_canonical_identifier() == other
        return (
            type(self) == type(other) and hash(self), hash(other)
        )

    def __str__(self):
        """Returns the SMILES string of the molecule"""
        return self.smiles

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class RingStack(object):
    """A class for holding ring information of a scaffold.
    This class is initialized by the Scaffold class.
    """

    __slots__ = (
        'owner',
        'info',
        'atom_rings',
        'bond_rings',
    )

    def __init__(self, owner):
        self.owner = weakref.proxy(owner)
        self.info = self.owner.mol.GetRingInfo()
        self.atom_rings = self.info.AtomRings()
        self.bond_rings = self.info.BondRings()

    @property
    def count(self):
        """Returns the number of rings in the stack"""
        return len(self)

    def __getitem__(self, index):
        return Ring(
            self.owner,
            self.atom_rings[index],
            self.bond_rings[index]
        )

    def __len__(self):
        """Returns the number of rings in the stack"""
        return self.info.NumRings()

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class Ring(object):
    """A class for holding information about a single ring.
    This class is initialized by the RingStack class.
    """

    __slots__ = 'owner', 'aix', 'bix'

    def __init__(self, owner, atom_indexes, bond_indexes):
        self.owner = owner
        self.aix = atom_indexes
        self.bix = bond_indexes

    @property
    def atoms(self):
        """Returns the rdkit atoms in this ring"""
        return [self.owner.atoms[x] for x in self.aix]

    @property
    def bonds(self):
        """Returns the rdkit bonds in this ring"""
        return [self.owner.bonds[x] for x in self.bix]

    @property
    def size(self):
        """Returns the size of the ring (number of atoms)"""
        return len(self)

    def __len__(self):
        """Returns the size of the ring (number of atoms)"""
        return len(self.aix)

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class RingSystemStack(object):
    """A class for holding the ring system information of a scaffold.
    This class is initialized by the Scaffold class.
    """

    __slots__ = (
        'owner',
        'atom_rings',
        'bond_rings',
        'ring_indexes'
    )

    def __init__(self, owner):
        self.owner = weakref.proxy(owner)
        systems = self._find_fused_rings()
        self.ring_indexes, self.atom_rings, self.bond_rings = systems

    def _find_fused_rings(self):
        """Private: find fused rings in the molecule"""
        rings = zip(self.owner.rings.atom_rings, self.owner.rings.bond_rings)
        rings = [([i], set(a), set(b)) for i, (a, b) in enumerate(rings)]
        systems = []
        for ring in rings:
            atoms = ring[1]
            n_systems = []
            for system in systems:
                common = atoms.intersection(system[1])
                if common:
                    ring[1].update(system[1])
                    ring[2].update(system[2])
                    ring[0].extend(system[0])
                else:
                    n_systems.append(system)
            n_systems.append(ring)
            systems = n_systems
        if systems:
            idx, aix, bix = zip(*systems)
            return tuple(map(tuple, idx)), tuple(map(tuple, aix)), tuple(map(tuple, bix))
        else:
            return (), (), ()

    @property
    def count(self):
        """Returns the number of ring systems in the stack"""
        return len(self)

    def __getitem__(self, index):
        return RingSystem(
            self.owner,
            self.atom_rings[index],
            self.bond_rings[index],
            self.ring_indexes[index]
        )

    def __len__(self):
        """Returns the number of ring systems in the stack"""
        return len(self.atom_rings)

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class RingSystem(object):
    """A class for holding information about a ring system.
    This class is initialized by the RingSystemStack class.
    """

    __slots__ = 'owner', 'aix', 'bix', 'rix'

    def __init__(self, owner, atom_indexes, bond_indexes, ring_indexes):
        self.owner = owner
        self.aix = atom_indexes
        self.bix = bond_indexes
        self.rix = ring_indexes

    @property
    def atoms(self):
        """Returns the rdkit atoms in this ring system"""
        return [self.owner.atoms[x] for x in self.aix]

    @property
    def bonds(self):
        """Returns the rdkit bonds in this ring system"""
        return [self.owner.bonds[x] for x in self.bix]

    @property
    def size(self):
        """Returns the size of the ring system (number of atoms)"""
        return len(self)

    def get_rings(self):
        """Return Ring objects which are subsets of the fused system"""
        return iter(self)

    def __getitem__(self, index):
        return self.owner.rings[self.rix[index]]

    def __len__(self):
        return len(self.aix)

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )
