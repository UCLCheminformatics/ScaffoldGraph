"""
scaffoldgraph.core.scaffold

A module defining scaffold objects used in ScaffoldGraph
"""

import weakref

from rdkit.Chem import MolToSmiles, BondType


class Scaffold(object):
    """A convenience class for defining scaffolds in scaffoldgraph.

    The class is fundamentally a wrapper around an rdkit Mol,
    containing methods and properties useful during fragmentation
    and prioritization. It is designed to be used internally by
    scaffoldgraph.

    Attributes
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The rdkit molecule representing the scaffold.
    hash_func : callable, None
        The hash function used to generate a canonical representation of
        the scaffold. If the `hash_func` is None (not set) then the
        canonical SMILES is used as the representation.

    Examples
    --------
    Initializing a Scaffold:

    >>> from rdkit import Chem
    >>> smiles = 'O=c1c(-c2ccccc2)cccn1Cc1c[nH]cn1'
    >>> mol = Chem.MolFromSmiles(smiles)
    >>> scaffold = Scaffold(mol)

    The underlying Mol can still be accessed:

    >>> scaffold.mol
    <rdkit.Chem.rdchem.Mol object at 0x7f03e94845d0>

    Accessing bonds & atoms is easy:

    >>> scaffold.atoms
    [<rdkit.Chem.rdchem.Atom at 0x7f03e9495940>,
    <rdkit.Chem.rdchem.Atom at 0x7f03e9495b70>,
    ...]
    >>> scaffold.bonds
    [<rdkit.Chem.rdchem.Bond at 0x7f03e9498170>,
    <rdkit.Chem.rdchem.Bond at 0x7f03e94982b0>,
    ...]

    Rings can also be accessed:

    >>> scaffold.rings
    <RingStack at 0x7f03e947cc60>
    >>> scaffold.rings.count
    3
    >>> ring = scaffold.rings[0]
    <Ring at 0x7f03e946ea08>
    >>> ring.size
    5
    >>> ring.aromatic
    True
    >>> atoms = ring.atoms
    >>> bonds = ring.bonds

    Scaffolds also holds information about ring systems and
    their properties:

    >>> scaffold.ring_systems
    <RingSystemStack at 0x7f03e948c750>
    >>> system = scaffold.ring_systems[0]
    >>> size = system.size
    >>> atoms = system.size
    >>> bonds = system.bonds

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
        """Initialize a Scaffold object.

        Parameters
        ----------
        rdmol: rdkit.Chem.rdchem.Mol
            An rdkit molecule representing a molecular scaffold. The
            class does not check the validity of the scaffold, as this
            definition may differ as per application, thus these checks
            should be performed by the user/application.
        hash_func: callable, optional
            Function for hash calculation, by default the canonical
            smiles string is used. The default is None. Ideally the
            hash provided should be canonical so that scaffolds which
            are identical according to a paticular definition are mapped
            to an identical hash. A user may want to customize the hash
            function to suit a specific equivalence definition.

        """
        self.mol = rdmol
        self._atoms = None
        self._bonds = None
        self._rings = None
        self._ring_systems = None
        self.hash_func = hash_func

    @property
    def name(self):
        """str : Returns the name/identifier of the scaffold"""
        try:
            return self.mol.GetProp('_Name')
        except KeyError:
            return None

    @name.setter
    def name(self, value):
        """Set the name/identifier of the scaffold."""
        self.mol.SetProp('_Name', str(value))

    @property
    def atoms(self):
        """list : Return a list of atoms in the scaffold."""
        if self._atoms is None:
            self._atoms = [a for a in self.mol.GetAtoms()]
        return self._atoms

    @property
    def bonds(self):
        """list : Return a list of bonds in the scaffold."""
        if self._bonds is None:
            self._bonds = [b for b in self.mol.GetBonds()]
        return self._bonds

    @property
    def rings(self):
        """RingStack : Return ring information for rings in the scaffold."""
        if self._rings is None:
            self._rings = RingStack(self)
        return self._rings

    @property
    def ring_systems(self):
        """RingSystemStack : Return ring system information for scaffold."""
        if self._ring_systems is None:
            self._ring_systems = RingSystemStack(self)
        return self._ring_systems

    @property
    def smiles(self):
        """str : Returns the canonical smiles string of the scaffold."""
        if not hasattr(self, '_smiles'):
            setattr(self, '_smiles', MolToSmiles(self.mol))
        return self._smiles

    @property
    def prioritization_rule(self):
        """str : Returns the prioritization rule used to select this scaffold."""
        try:
            return self.mol.GetProp('prioritization_rule')
        except KeyError:
            return None

    @prioritization_rule.setter
    def prioritization_rule(self, value):
        """Sets the prioritization rule used to select this scaffold."""
        self.mol.SetProp('prioritization_rule', str(value))

    @property
    def removed_ring_idx(self):
        """int : Returns the index of the removed ring in the child molecule."""
        try:
            return self.mol.GetIntProp('removed_ring_idx')
        except KeyError:
            return None

    @removed_ring_idx.setter
    def removed_ring_idx(self, value):
        """Sets the index of the removed ring in the child molecule."""
        self.mol.SetIntProp('removed_ring_idx', int(value))

    def get_canonical_identifier(self):
        """Returns a canonical identifier for the scaffold.

        The canonical identifier is determined by the internal
        `hash_func` attribute. This can be set upon initialization
        or afterwards by: Scaffold.hash_func = callable.

        Returns
        -------
        str
            A canonical identifier for the scaffold.

        """
        if self.hash_func:
            return self.hash_func(self.mol)
        return self.smiles

    def __bool__(self):
        """Returns True if the molecule contains at least 1 atom."""
        return len(self.atoms) >= 1

    def __hash__(self):
        """Returns a hash of the canonical identifier."""
        return hash(self.get_canonical_identifier())

    def __eq__(self, other):
        if isinstance(other, str):
            return self.get_canonical_identifier() == other
        return (
            type(self) == type(other) and hash(self), hash(other)
        )

    def __str__(self):
        """Returns the SMILES string of the molecule."""
        return self.smiles

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class RingStack(object):
    """A class for holding ring information of a Scaffold.

    This class is initialized by the Scaffold class.

    Attributes
    ----------
    owner : weakproxy (Scaffold)
        A weak reference to the owning Scaffold object.
    info : rdkit.Chem.rdchem.RingInfo
        Ring information for the scaffolds rings.
    atom_rings : tuple
        A tuple of tuples containing the atom indicies for
        each ring.
    bond_rings : tuple
        A tuple of tuples containing the bond indicies for
        each ring.

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
        """int : Returns the number of rings in the stack."""
        return len(self)

    def __getitem__(self, index):
        return Ring(
            self.owner,
            self.atom_rings[index],
            self.bond_rings[index]
        )

    def __len__(self):
        """Returns the number of rings in the stack."""
        return self.info.NumRings()

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class Ring(object):
    """A class for holding information about a single ring.

    This class is initialized by the RingStack class.

    Attributes
    ----------
    owner : weakproxy (Scaffold)
        A weak reference to the owning Scaffold object.
    aix : tuple
        Indicies of atoms in the ring.
    bix : tuple
        Indicies of bonds in the ring.

    """
    __slots__ = 'owner', 'aix', 'bix'

    def __init__(self, owner, atom_indexes, bond_indexes):
        self.owner = owner
        self.aix = atom_indexes
        self.bix = bond_indexes

    @property
    def atoms(self):
        """list : Returns the rdkit atoms in this ring."""
        return [self.owner.atoms[x] for x in self.aix]

    @property
    def bonds(self):
        """list : Returns the rdkit bonds in this ring."""
        return [self.owner.bonds[x] for x in self.bix]

    @property
    def size(self):
        """int : Returns the size of the ring (number of atoms)."""
        return len(self)

    @property
    def aromatic(self):
        """bool : Return True/False if ring is/isn't aromatic."""
        return all([x.GetIsAromatic() for x in self.bonds])

    def get_attachment_points(self, include_exocyclic=False):
        """Return atom indexes for linker attachments.

        Parameters
        ----------
        include_exocyclic : bool, optional
            Whether to include exocyclic double-bonded
            attachments. The default is False.

        Returns
        -------
        set
            Atom indicies, of atoms which are attachment
            points for linkers.

        """
        attachments = set()
        ri = self.owner.rings.info
        for atom in self.atoms:
            index = atom.GetIdx()
            if (
                    ri.NumAtomRings(index) == 1 or
                    any([not b.IsInRing() for b in atom.GetBonds()])
            ):
                if atom.GetDegree() > 2:
                    if (
                        include_exocyclic is False
                        and self.is_exocyclic_attachment(atom)
                    ):
                        continue
                    attachments.add(index)
        return attachments

    def is_exocyclic_attachment(self, atom):
        """Returns if a ring atom is an exocyclic attachment point.

        Parameters
        ----------
        atom : rdkit.Chem.rdchem.Atom
            The ring atom to test for exocyclic attachment.

        Returns
        -------
        bool
            Returns True if the supplied atom is attached to
            an exocyclic attachment.

        Raises
        ------
        ValueError
            If the supplied atom is not within the ring.

        """
        exocyclic = False
        if atom not in self.atoms:
            raise ValueError(f'atom {atom.GetIdx()} not in ring')
        ri = self.owner.rings.info
        for bond in atom.GetBonds():
            if ri.NumBondRings(bond.GetIdx()) == 0:
                other = bond.GetOtherAtom(atom)
                if (
                    bond.GetBondType() == BondType.DOUBLE
                    and other.GetDegree() == 1
                ):
                    exocyclic = True
        return exocyclic

    def get_ring_system(self):
        """Return the ring system associated with this ring.

        Returns
        -------
        RingSystem
            The ring system containing this ring. Note that
            the returned ring system may only contain this
            ring if the system is of `num_rings` == 1.

        See Also
        --------
        RingSystem

        """
        for system in self.owner.ring_systems:
            if system.contains(self):
                return system

    def __len__(self):
        """Returns the size of the ring (number of atoms)."""
        return len(self.aix)

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class RingSystemStack(object):
    """A class for holding the ring system information of a Scaffold.

    This class is initialized by the Scaffold class.

    Attributes
    ----------
    owner : weakproxy (Scaffold)
        A weak reference to the owning Scaffold object.
    ring_indexes : tuple
        A tuple of tuples containg the ring indicies for
        each ring system.
    atom_rings : tuple
        A tuple of tuples containing the atom indicies for
        each ring system.
    bond_rings : tuple
        A tuple of tuples containing the bond indicies for
        each ring system.

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
        """Private: find fused rings in the molecule."""
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
        """int : Returns the number of ring systems in the stack."""
        return len(self)

    def __getitem__(self, index):
        return RingSystem(
            self.owner,
            self.atom_rings[index],
            self.bond_rings[index],
            self.ring_indexes[index]
        )

    def __len__(self):
        """Returns the number of ring systems in the stack."""
        return len(self.atom_rings)

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )


class RingSystem(object):
    """A class for holding information about a ring system.

    This class is initialized by the RingSystemStack class.

    Attributes
    ----------
    owner : weakproxy (Scaffold)
        A weak reference to the owning Scaffold object.
    aix : tuple
        Indicies of atoms in the ring system.
    bix : tuple
        Indicies of bonds in the ring sytem.
    rix : tuple
        Indicies of rings in the ring system.

    """
    __slots__ = 'owner', 'aix', 'bix', 'rix'

    def __init__(self, owner, atom_indexes, bond_indexes, ring_indexes):
        self.owner = owner
        self.aix = atom_indexes
        self.bix = bond_indexes
        self.rix = ring_indexes

    @property
    def atoms(self):
        """list : Returns the rdkit atoms in this ring system."""
        return [self.owner.atoms[x] for x in self.aix]

    @property
    def bonds(self):
        """list : Returns the rdkit bonds in this ring system."""
        return [self.owner.bonds[x] for x in self.bix]

    @property
    def size(self):
        """int : Returns the size of the ring system (number of atoms)."""
        return len(self)

    @property
    def num_rings(self):
        """int : Return the number of rings in the ring system."""
        return len(self.rix)

    def get_rings(self):
        """Return Ring objects which are subsets of the ring system.

        Returns
        -------
        iterable
            An iterable of Ring objects which are subsets of this
            ring system.
        """
        return iter(self)

    def contains(self, ring):
        """Return True if ring system contains query ring.

        Parameters
        ----------
        ring : Ring
            A ring query.

        Returns
        -------
        bool
            True if the query ring is contained in this ring
            system.

        """
        assert type(ring) == Ring, f'query must be a {Ring} object'
        return len(set(self.aix).intersection(ring.aix)) > 0

    def contains_ring_idx(self, ring_idx):
        """Return True if ring system contains query ring index.

        Parameters
        ----------
        ring_idx : int
            Index of query ring.

        Returns
        -------
        bool
            True if ring with the supplied index is contained
            within this ring system.

        """
        return ring_idx in self.rix

    def get_attachment_points(self, include_exocyclic=False):
        """Return atom indexes for linker attachments.

        Parameters
        ----------
        include_exocyclic : bool, optional
            Whether to include exocyclic double-bonded
            attachments. The default is False.

        Returns
        -------
        set
            Atom indicies, of atoms which are attachment
            points for linkers.

        """
        attachments = set()
        for atom in self.atoms:
            index = atom.GetIdx()
            if any([not b.IsInRing() for b in atom.GetBonds()]):
                if atom.GetDegree() > 2:
                    if (
                        include_exocyclic is False
                        and self.is_exocyclic_attachment(atom)
                    ):
                        continue
                    attachments.add(index)
        return attachments

    def is_exocyclic_attachment(self, atom):
        """Returns if a ring atom is an exocyclic attachment point.

        Parameters
        ----------
        atom : rdkit.Chem.rdchem.Atom
            The ring atom to test for exocyclic attachment.

        Returns
        -------
        bool
            Returns True if the supplied atom is attached to
            an exocyclic attachment.

        Raises
        ------
        ValueError
            If the supplied atom is not within the ring system.

        """
        exocyclic = False
        if atom not in self.atoms:
            raise ValueError(f'atom {atom.GetIdx()} not in ring system')
        ri = self.owner.rings.info
        for bond in atom.GetBonds():
            if ri.NumBondRings(bond.GetIdx()) == 0:
                other = bond.GetOtherAtom(atom)
                if (
                    bond.GetBondType() == BondType.DOUBLE
                    and other.GetDegree() == 1
                ):
                    exocyclic = True
        return exocyclic

    def __getitem__(self, index):
        return self.owner.rings[self.rix[index]]

    def __len__(self):
        return len(self.aix)

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )
