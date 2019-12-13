"""
scaffoldgraph.core.graph

Defines the base class for scaffold graphs in ScaffoldGraph
"""

from abc import ABC, abstractmethod
from collections import Counter

import tqdm
from loguru import logger
from networkx import DiGraph
from rdkit import RDLogger
from rdkit.Chem import rdMolHash, MolToSmiles, rdmolops
from rdkit.Chem.rdMolDescriptors import CalcNumRings

from scaffoldgraph.io import *
from .fragment import get_murcko_scaffold, get_annotated_murcko_scaffold
from .scaffold import Scaffold

rdlogger = RDLogger.logger()


def init_molecule_name(mol):
    """Initialize the name of a molecule if not provided"""
    if not mol.HasProp('_Name') or mol.GetProp('_Name') == '':
        n = rdMolHash.GenerateMoleculeHashString(mol)
        mol.SetProp('_Name', n)


class ScaffoldGraph(DiGraph, ABC):
    """Abstract base class for ScaffoldGraphs"""

    def __init__(self, graph=None, fragmenter=None):
        """
        Initialize a ScaffoldGraph object

        Parameters
        ----------
        graph: Graph data to inherit (optional, default=None)
        """
        super(ScaffoldGraph, self).__init__(graph)
        self.fragmenter = fragmenter

    def _construct(self, molecules, ring_cutoff=10, progress=False, annotate=True):
        """
        Private method for graph construction, called by constructors

        Parameters
        ----------
        molecules: iterable of rdkit molecules for processing
        ring_cutoff: ignore molecules with more than the specified number
            of rings to avoid extended processing times (default: 10)
        annotate: if True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold)
        progress: if True show a progress bar monitoring progress (default: False)
        """
        rdlogger.setLevel(4)  # Suppress the RDKit logs
        progress = progress is False
        desc = self.__class__.__name__
        for molecule in tqdm.tqdm(molecules, disable=progress, desc=desc, miniters=1, dynamic_ncols=True):
            if molecule is None:  # logged in suppliers
                continue
            init_molecule_name(molecule)
            if CalcNumRings(molecule) > ring_cutoff:
                name = molecule.GetProp('_Name')
                logger.warning(f'Molecule {name} filtered (> {ring_cutoff} rings)')
                continue
            rdmolops.RemoveStereochemistry(molecule)
            scaffold = Scaffold(get_murcko_scaffold(molecule))
            if scaffold:  # Checks that a scaffold has at least 1 atom
                if annotate:
                    annotation = get_annotated_murcko_scaffold(molecule, scaffold.mol, False)
                else:
                    annotation = None
                self.add_scaffold_node(scaffold)
                self.add_molecule_node(molecule)
                self.add_molecule_edge(molecule, scaffold, annotation=annotation)
                if scaffold.rings.count > 1:
                    self._recursive_constructor(scaffold)
            else:
                name = molecule.GetProp('_Name')
                logger.warning(f'No top level scaffold for molecule {name}')
        rdlogger.setLevel(3)  # Enable the RDKit logs

    @abstractmethod
    def _recursive_constructor(self, child):
        """This method should be implemented by the subclass, used during recursive
        fragmentation of a scaffold"""
        raise NotImplementedError()

    @property
    def num_scaffold_nodes(self):
        """Return the number of scaffolds in the graph"""
        return len(list(self.get_scaffold_nodes()))

    @property
    def num_molecule_nodes(self):
        """Return the number of molecules in the graph"""
        return len(list(self.get_molecule_nodes()))

    def get_scaffold_nodes(self, data=False):
        if data is True:
            return ((n, self.nodes[n]) for n, d in self.nodes(data='type') if d == 'scaffold')
        else:
            return (n for n, d in self.nodes(data='type') if d == 'scaffold')

    def get_molecule_nodes(self, data=False):
        if data is True:
            return ((n, self.nodes[n]) for n, d in self.nodes(data='type') if d == 'molecule')
        else:
            return (n for n, d in self.nodes(data='type') if d == 'molecule')

    def get_hierarchy_sizes(self):
        hierarchy = (d['hierarchy'] for _, d in self.get_scaffold_nodes(data=True))
        return Counter(hierarchy)

    def add_molecule_node(self, molecule, **attr):
        name = molecule.GetProp('_Name')
        default_attr = dict(type='molecule', smiles=MolToSmiles(molecule))
        default_attr.update(attr)
        self.add_node(name, **default_attr)

    def add_scaffold_node(self, scaffold, **attr):
        default_attr = dict(type='scaffold', hierarchy=scaffold.rings.count)
        default_attr.update(attr)
        self.add_node(scaffold.get_canonical_identifier(), **default_attr)

    def add_molecule_edge(self, molecule, scaffold, **attr):
        name = molecule.GetProp('_Name')
        default_attr = dict(type=0)
        default_attr.update(attr)
        self.add_edge(scaffold.get_canonical_identifier(), name, **default_attr)

    def add_scaffold_edge(self, parent, child, **attr):
        default_attr = dict(type=1)
        default_attr.update(attr)
        self.add_edge(parent.get_canonical_identifier(),
                      child.get_canonical_identifier(),
                      **default_attr)

    @classmethod
    def from_sdf(cls, file_name, ring_cutoff=10, progress=False, annotate=True, **kwargs):
        """Construct a scaffoldgraph from an SDF file.

        Parameters
        ----------
        file_name (str): path to an SDF input file
        ring_cutoff (int): ignore molecules with more rings than this cutoff (default: 10)
        progress (bool): if True show a progress bar to monitor construction progress (default: False)
        annotate: if True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold)
        kwargs: arguments to pass to the scaffoldgraph __init__
        """

        with open(file_name, 'rb') as sdf:
            supplier = read_sdf(sdf, requires_length=progress is True)
            instance = cls(**kwargs)
            instance._construct(supplier, ring_cutoff=ring_cutoff, progress=progress, annotate=annotate)
        return instance

    @classmethod
    def from_smiles_file(cls, file_name, delimiter=' ', smiles_column=0, name_column=1,
                         header=False, ring_cutoff=10, progress=False, annotate=True, **kwargs):
        """Construct a scaffoldgraph from a SMILES file.

        Parameters
        ----------
        file_name (str): path to a SMILES input file
        delimiter (str): delimiter used in SMILES file (default: ' ')
        smiles_column (int): index of column containing SMILES strings (default: 0)
        name_column (int): index of column containing molecule names (default: 1)
        header (bool): if True skip the first line of the SMILES file
        ring_cutoff (int): ignore molecules with more rings than this cutoff (default: 10)
        progress (bool): if True show a progress bar to monitor construction progress (default: False)
        annotate: if True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold)
        kwargs: arguments to pass to the scaffoldgraph __init__
        """

        supplier = read_smiles_file(file_name, delimiter, smiles_column, name_column,
                                    header, requires_length=progress is True)
        instance = cls(**kwargs)
        instance._construct(supplier, ring_cutoff=ring_cutoff, progress=progress, annotate=annotate)
        return instance

    @classmethod
    def from_supplier(cls, supplier, ring_cutoff=10, progress=False, annotate=True, **kwargs):
        """Construct a scaffoldgraph from a custom Mol supplier

        Parameters
        ----------
        supplier (iterable): A custom supplier of rdkit molecules, it is expected that each molecule
            will be assigned a property '_Name' serving as an identifier for that molecule. To use
            the progress option the supplier must have an associated length.
        ring_cutoff (int): ignore molecules with more rings than this cutoff (default: 10)
        progress (bool): if True show a progress bar to monitor construction progress (default: False)
        annotate: if True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold)
        kwargs: arguments to pass to the scaffoldgraph __init__
        """

        instance = cls(**kwargs)
        instance._construct(supplier, ring_cutoff=ring_cutoff, progress=progress, annotate=annotate)
        return instance

    @classmethod
    def from_dataframe(cls, df, smiles_column='Smiles', name_column='Name', ring_cutoff=10,
                       progress=False, annotate=True, **kwargs):
        """Construct a scaffoldgraph from a pandas DataFrame.

        Parameters
        ----------
        df (pd.DataFrame): a pandas DataFrame containing SMILES strings and a molecule identifier
        smiles_column (str): label of column containing SMILES strings (default: 'Smiles')
        name_column (str): label of column containing SMILES strings (default: 'Name')
        ring_cutoff (int): ignore molecules with more rings than this cutoff (default: 10)
        progress (bool): if True show a progress bar to monitor construction progress (default: False)
        annotate: if True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold)
        """

        supplier = read_dataframe(df, smiles_column, name_column)
        instance = cls(**kwargs)
        instance._construct(supplier, ring_cutoff=ring_cutoff, progress=progress, annotate=annotate)
        return instance

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )
