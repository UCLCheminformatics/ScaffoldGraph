"""
scaffoldgraph.core.graph

Defines the base class for scaffold graphs in ScaffoldGraph
"""

from abc import ABC, abstractmethod
from collections import Counter

import networkx as nx
import gzip

from loguru import logger
from tqdm.auto import tqdm

from rdkit import __version__ as rdversion
from rdkit.Chem import rdMolHash, MolToSmiles, rdmolops
from rdkit.Chem.rdMolDescriptors import CalcNumRings

from scaffoldgraph.io import *
from scaffoldgraph.utils import canonize_smiles, suppress_rdlogger

from .fragment import (
    get_murcko_scaffold,
    get_annotated_murcko_scaffold,
    keep_largest_fragment,
    flatten_isotopes,
    discharge_and_deradicalize,
)

from .scaffold import Scaffold


def init_molecule_name(mol):
    """Initialize the name of a molecule if not provided.

    If the molecule has no `_Name` property then it is
    set as a hash computed by ``rdkit.Chem.rdMolHash``.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol

    Notes
    -----
    Since rdkit 2020.09.01 the GetMolHashString has been
    deprecated. If using an rdkit version >= 2020.09.01
    the function sets the name to 'MolNode-' plus the
    canonical smiles. This prevents collisions with
    scaffolds which are hashed using their canonical smiles.

    """
    if not mol.HasProp('_Name') or mol.GetProp('_Name') == '':
        if rdversion < '2020.09.01':
            n = rdMolHash.GenerateMoleculeHashString(mol)
        else:  # New version deprecated GenrateMolHashString
            hashf = rdMolHash.HashFunction.CanonicalSmiles
            n = 'MolNode-' + rdMolHash.MolHash(mol, hashf)
        mol.SetProp('_Name', n)


class ScaffoldGraph(nx.DiGraph, ABC):
    """Base class for ScaffoldGraphs.

    A ScaffoldGraph stores nodes and edges with optional data, or attributes.
    Nodes may represent scaffolds or molecules (although other node types may
    also be added). Scaffold nodes are keyed by a hash, this is usually a
    canonical SMILES string. To avoid molecule node keys confiliciting with
    scaffold node keys (share the same canonical SMILES), they are keyed by a
    name provided by the molecule property '_Name'. If this property is not
    available a hash is generated for the molecule which differs from the
    canonical SMILES.

    ScaffoldGraphs hold directed edges. An edge represents a hierarchial
    relationship between a molecule and a scaffold or a scaffold and another
    (other edge types may also exist). Edges may also contain data or attributes.

    Attributes
    ----------
    fragmenter : scaffoldgraph.core.fragment.Fragmenter
        A ``scaffoldgraph.core.fragment.Fragmenter`` class for producing
        the next scaffold set for a given molecular input.

    **Subclasses:**

    In scaffoldgraph all ScaffoldGraphs (e.g. ``ScaffoldNetwork``) inherit from
    the ``ScaffoldGraph`` base class. It provides an interface to simplify
    building and querying hierarchial scaffold relationships. Each subclass
    represents a different way of constructing the graph. This can be customized
    in a few ways:

        Implementing the _hierarchy_constructor function. This function
        is responsible for taking the next set of fragments and building
        the graph structure (subclasses must implment this function). For
        example in the ScaffoldTree class the funtion uses a prioritization
        system for determining which scaffold to add to the graph. This
        function MUST be implemented.

        By providing a ``Fragmenter`` callable to produce the next set of
        scaffolds to be added to the graph, the way in which scaffolds are
        fragmented can be customized.

        By mofifying scaffold initilaization. Scaffold initilization can be customized
        by mofifying the functions; _initialize_scaffold and _preprocess_scaffold.
        _initialize scaffold is responsible for creating the top-level scaffold
        and adding the molecule node, scaffold node and connecting these via an edge.
        The _initialize_scaffold function also calls the _preprocess_scaffold function
        which contains optional methods for preprocessing scaffolds, i.e. removing
        specific isotopes.

    See Also
    --------
    networkx.classes.digraph.DiGraph
    ScaffoldNetwork
    ScaffoldTree
    HierS

    """
    def __init__(self, graph=None, fragmenter=None, graph_type=None, **attr):
        """ Initialize a ScaffoldGraph.

        Parameters
        ----------
        graph : input graph, optional
            Data to initialize graph. If None (default) an empty
            graph is created. The data can be any format that is supported
            by the ``to_networkx_graph()`` function, currently including
            edge list, dict of dicts, dict of lists, NetworkX graph,
            NumPy matrix or 2d ndarray, SciPy sparse matrix,
            or PyGraphviz graph. This argument is passed to the networkx
            DiGraph constructor.
        fragmenter : scaffoldgraph.core.fragment.Fragmenter
            A ``scaffoldgraph.core.fragment.Fragmenter`` class for producing
            the next scaffold set for a given molecular input.
        graph_type : str, optional
            The type of graph being constructed (set as a global
            graph attribute).
        **attr : keyword arguments, optional
            Attributes to add to graph as key=value pairs. The default is
            no attributes.

        """
        super(ScaffoldGraph, self).__init__(graph, graph_type=graph_type, **attr)
        self.graph.setdefault('num_linear', 0)    # track number of linear mols.
        self.graph.setdefault('num_filtered', 0)  # track number of filtered mols.
        self.fragmenter = fragmenter

    @suppress_rdlogger()
    def _construct(self, molecules, init_args, ring_cutoff=10, progress=False):
        """Private method for graph construction, called by constructors.

        The constructor is fairly generic allowing the user to customise
        hierarchy construction by changing the _hierarchy_constructor
        function. The user also has further control of this process and
        is able to change how a scaffold is initialized (`_initialize_scaffold`),
        how it is preprocessed (`_preprocess_scaffold`) and how molecules with
        no top-level scaffold (i.e. linear molecules) are handled.

        Parameters
        ----------
        molecules : iterable
            An iterable of rdkit molecules for processing
        init_args : dict
            A dictionary containing arguments for scaffold initialization and
            preprocessing.
        ring_cutoff : int, optional
            Ignore molecules with more than the specified number of rings to avoid
            extended processing times. The default is 10.
        progress : bool, optional
            If True show a progress bar monitoring progress. The default is False.

        See Also
        --------
        _initialize_scaffold
        _preprocess_scaffold
        _process_no_top_level

        """
        desc, progress = self.__class__.__name__, progress is False
        for molecule in tqdm(
            molecules, disable=progress, desc=desc,
            miniters=1, dynamic_ncols=True,
        ):
            if molecule is None:
                continue
            init_molecule_name(molecule)
            if CalcNumRings(molecule) > ring_cutoff:
                name = molecule.GetProp('_Name')
                logger.warning(f'Molecule {name} filtered (> {ring_cutoff} rings)')
                self.graph['num_filtered'] = self.graph.get('num_filtered', 0) + 1
                continue
            scaffold = self._initialize_scaffold(molecule, init_args)
            if scaffold is not None:
                self._hierarchy_constructor(scaffold)

    def _initialize_scaffold(self, molecule, init_args):
        """Initialize the top-level scaffold for a molecule.

        Initialization generates a murcko scaffold, performs
        any preprocessing required and then adds the scaffold
        node to the graph connecting it to its child molecule.
        This process can be customised in subclasses to modify
        how a scaffold is initialized.

        Parameters
        ----------
        molecule : rdkit.Chem.rdchem.Mol
            A molecule from whicg to initialize a scaffold.
        init_args : dict
            A dictionary containing arguments for scaffold
            initialization and preprocessing.

        Returns
        -------
        scaffold : Scaffold
            A Scaffold object containing the initialized
            scaffold to be processed further (hierarchy
            generation).

        """
        scaffold_rdmol = get_murcko_scaffold(molecule)
        if scaffold_rdmol.GetNumAtoms() <= 0:
            return self._process_no_top_level(molecule)
        scaffold_rdmol = self._preprocess_scaffold(scaffold_rdmol, init_args)
        scaffold = Scaffold(scaffold_rdmol)
        annotation = None
        if init_args.get('annotate') is True:
            annotation = get_annotated_murcko_scaffold(molecule, scaffold_rdmol, False)
        self.add_scaffold_node(scaffold)
        self.add_molecule_node(molecule)
        self.add_molecule_edge(molecule, scaffold, annotation=annotation)
        return scaffold

    def _preprocess_scaffold(self, scaffold_rdmol, init_args):
        """Preprocess a scaffold before initialization.

        When subclassing a user is able to customize/extend
        preprocessing options to suit their needs. The default
        method has the option to; remove stereochemistry,
        remove specific isotopes, keep the largest disconnected
        fragment and discharge and deradicalize a molecule. The
        options are controlled by the init_args dictionary. This
        argument dictionary is formed by the constructors. Note
        that the stereochemistry argument is specific to the
        hierarchy generating method so isn't specified in the
        constructors.

        Parameters
        ----------
        scaffold_rdmol : rdkit.Chem.rdchem.Mol
            rdkit molecule of scaffold to be preprocessed
        init_args : dict
            A dictionary containing arguments for scaffold initialization
            and preprocessing.

        Returns
        -------
        scaffold_rdmol : rdkit.Chem.rdchem.Mol
            rdkit molecule containing preprocessed scaffold.

        """
        if init_args.get('remove_stereo', True) is True:
            rdmolops.RemoveStereochemistry(scaffold_rdmol)
        if init_args.get('flatten_isotopes', False) is True:
            flatten_isotopes(scaffold_rdmol)
        if init_args.get('keep_largest', False) is True:
            scaffold_rdmol = keep_largest_fragment(scaffold_rdmol)
        if init_args.get('discharge', False) is True:
            scaffold_rdmol = discharge_and_deradicalize(scaffold_rdmol)
        return scaffold_rdmol

    def _process_no_top_level(self, molecule):
        """Private: Process molecules with no top-level scaffold.

        Can be overwritten to change the behaviour when a molecule
        does not contain a top-level scaffold (i.e. is linear). By
        default a warning is logged and the molecule is counted to
        the num_linear graph attribute (`self.graph['num_linear']`).
        The molecule node is not added to the graph.

        If for example a user wanted to track what molecules do not
        have a top level scaffold, they could overide this function
        to append the names of these molecules to a list, or indeed
        add them to the graph.

        Parameters
        ----------
        molecule : rdkit.Chem.rdchem.Mol
            An rdkit molecule determined to have no top-level scaffold.

        """
        name = molecule.GetProp('_Name')
        logger.warning(f'No top level scaffold for molecule: {name}')
        self.graph['num_linear'] = self.graph.get('num_linear', 0) + 1
        return None

    @abstractmethod
    def _hierarchy_constructor(self, child):
        """
        This method should be implemented by the subclass, used during
        fragmentation of a scaffold to generate a scaffold hierarchy.

        Parameters
        ----------
        child : scaffoldgraph.core.Scaffold

        """
        raise NotImplementedError()

    @property
    def num_scaffold_nodes(self):
        """int : Return the number of scaffold nodes in the graph."""
        count = 0
        for _ in self.get_scaffold_nodes():
            count += 1
        return count

    @property
    def num_molecule_nodes(self):
        """int : Return the number of molecule nodes in the graph."""
        count = 0
        for _ in self.get_molecule_nodes():
            count += 1
        return count

    def _get_nodes_with_type(self, _type, data, default):
        """
        Private: Return a generator of all nodes in the graph with a 'type'
        attribute equal to _type.

        Parameters
        ----------
        _type : str
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.

        Returns
        -------
        nodes : generator
            A generator containing all nodes with a 'type' attribute equal
            to `_type`.

        """
        if data is False:
            return (n for n, d in self.nodes(data='type') if d == _type)
        elif data is True:
            return ((n, self.nodes[n]) for n, d in self.nodes(data='type') if d == _type)
        else:
            return ((n, self.nodes[n].get(data, default)) for n, d in self.nodes(data='type') if d == _type)

    def get_scaffold_nodes(self, data=False, default=None):
        """Return a generator of all scaffold nodes in the graph.

        Parameters
        ----------
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.

        Returns
        -------
        nodes : generator
            A generator containing all nodes with a 'type' attribute equal
            to `_type`.

        """
        return self._get_nodes_with_type('scaffold', data, default)

    def get_molecule_nodes(self, data=False, default=None):
        """Return a generator of all molecule nodes in the graph.

        Parameters
        ----------
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.

        Returns
        -------
        nodes : generator
            A generator containing all nodes with a 'type' attribute equal
            to `_type`.

        """
        return self._get_nodes_with_type('molecule', data, default)

    def get_hierarchy_sizes(self):
        """
        Return a ``collections.Counter`` object indicating the number of scaffolds
        within each hierarchy level.

        Returns
        -------
        collections.Counter
            A dict like object with keys representing hierarchies and values
            representing the number of scaffolds within the corresponding
            hierarchy.

        """
        hierarchy = (d['hierarchy'] for _, d in self.get_scaffold_nodes(data=True))
        return Counter(hierarchy)

    def max_hierarchy(self):
        """int : Return the largest hierarchy level"""
        return max(self.get_hierarchy_sizes())

    def min_hierarchy(self):
        """int : Return the smallest hierarchy level"""
        return min(self.get_hierarchy_sizes())

    def get_scaffolds_in_hierarchy(self, hierarchy):
        """Return a generator of all scaffolds within a specified hierarchy.

        Parameters
        ----------
        hierarchy : int
            The hierarchy level to retrieve.

        Returns
        -------
        nodes : generator
            A generator containing all  scaffold nodes belonging to the
            specified hierarchy.

        """
        for s, d in self.get_scaffold_nodes(data=True):
            if d['hierarchy'] == int(hierarchy):
                yield s

    def scaffold_in_graph(self, scaffold_smiles):
        """Returns True if the specified scaffold SMILES is in the scaffold graph.

        Parameters
        ----------
        scaffold_smiles : str
            SMILES of query scaffold.

        Returns
        -------
        bool
            True if the scaffold is in the graph.

        Notes
        -----
        If not initially found the SMILES is canonized and the graph is searched
        with the canonized SMILES key.

        """
        result = scaffold_smiles in self
        if result is not True:
            scaffold_smiles = canonize_smiles(scaffold_smiles, failsafe=True)
            result = scaffold_smiles in self
        return result and self.nodes[scaffold_smiles]['type'] == 'scaffold'

    def molecule_in_graph(self, molecule_id):
        """Returns True if specified molecule ID is in the scaffold graph.

        Parameters
        ----------
        molecule_id : str
            ID of query molecule.

        Returns
        -------
        bool
            True if the molecule is in the graph.

        """
        return str(molecule_id) in self and self.nodes[molecule_id]['type'] == 'molecule'

    def get_molecules_for_scaffold(self, scaffold_smiles, data=False, default=None):
        """Return a list of molecule IDs which are represented by a scaffold in the graph.

        Parameters
        ----------
        scaffold_smiles : str
            SMILES of query scaffold.
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.

        Returns
        -------
        list
            A list of molecule nodes.

        Notes
        -----
        Molecules are found by traversing the graph. In the case of a scaffold tree
        the results represent the rules used to prioritize the scaffolds.

        """
        molecules = []
        if scaffold_smiles not in self:
            scaffold_smiles = canonize_smiles(scaffold_smiles, failsafe=True)
            if scaffold_smiles not in self:
                return molecules
        for succ in nx.bfs_tree(self, scaffold_smiles, reverse=False):
            if self.nodes[succ].get('type') == 'molecule':
                if data is False:
                    molecules.append(succ)
                elif data is True:
                    molecules.append((succ, self.nodes[succ]))
                else:
                    molecules.append((succ, self.nodes[succ].get(data, default)))
        return molecules

    def get_scaffolds_for_molecule(self, molecule_id, data=False, default=None):
        """Return a list of scaffold SMILES connected to a query molecule ID.

        Parameters
        ----------
        molecule_id:  str
            ID of query molecule.
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.

        Returns
        -------
        list
            A list of scaffold nodes.

        Notes
        -----
        Scaffolds are found by traversing the graph. In the case of a scaffold tree
        the results represent the rules used to prioritize the scaffolds.

        """
        scaffolds = []
        if molecule_id not in self:
            return scaffolds
        for succ in nx.bfs_tree(self, molecule_id, reverse=True):
            if self.nodes[succ].get('type') == 'scaffold':
                if data is False:
                    scaffolds.append(succ)
                elif data is True:
                    scaffolds.append((succ, self.nodes[succ]))
                else:
                    scaffolds.append((succ, self.nodes[succ].get(data, default)))
        return scaffolds

    def _get_scaffold_hierarchy(self, scaffold_smiles, data=False, default=None, max_levels=-1, traversal='parent'):
        """Private: Return a list of parent/child scaffolds for a query scaffold.

        Parameters
        ----------
        scaffold_smiles : str
            SMILES of query scaffold.
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.
        max_levels : int, optional
            If > 0 only return scaffolds with a hierarchy difference to the
            query scaffold of `max_levels`.
        traversal : {'parent', 'child'}, optional
            Direction of traversal.

        Returns
        -------
        list
            A list of scaffold parent/child nodes.

        """
        assert traversal in {'parent', 'child'}
        reverse = traversal == 'parent'
        next_hiers = []
        if scaffold_smiles not in self:
            scaffold_smiles = canonize_smiles(scaffold_smiles, failsafe=True)
            if scaffold_smiles not in self:
                return next_hiers
        level = self.nodes[scaffold_smiles].get('hierarchy', float('inf'))
        bfs = iter(nx.bfs_tree(self, scaffold_smiles, reverse=reverse).nodes)
        next(bfs)  # first entry is the query node
        for succ in bfs:
            d = self.nodes[succ]
            if d.get('type') == 'scaffold' and (max_levels < 0 or level - d.get('hierarchy', 0) <= max_levels):
                if data is False:
                    next_hiers.append(succ)
                elif data is True:
                    next_hiers.append((succ, self.nodes[succ]))
                else:
                    next_hiers.append((succ, self.nodes[succ].get(data, default)))
        return next_hiers

    def get_parent_scaffolds(self, scaffold_smiles, data=False, default=None, max_levels=-1):
        """Return a list of parent scaffolds for a query scaffold.

        Parameters
        ----------
        scaffold_smiles : str
            SMILES of query scaffold.
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.
        max_levels : int, optional
            If > 0 only return scaffolds with a hierarchy difference to the
            query scaffold of `max_levels`.

        Returns
        -------
        list
            A list of scaffold parent nodes.

        """
        return self._get_scaffold_hierarchy(scaffold_smiles, data, default, max_levels, 'parent')

    def get_child_scaffolds(self, scaffold_smiles, data=False, default=None, max_levels=-1):
        """Return a list of child scaffolds for a query scaffold.

        Parameters
        ----------
        scaffold_smiles : str
            SMILES of query scaffold.
        data : str, bool, optional
            The scaffold node attribute returned in 2-tuple (n, ddict[data]).
            If True, return entire node attribute dict as (n, ddict).
            If False, return just the nodes n. The default is False.
        default : value, bool, optional
            Value used for nodes that don't have the requested attribute.
            Only relevant if data is not True or False.
        max_levels : int, optional
            If > 0 only return scaffolds with a hierarchy difference to the
            query scaffold of `max_levels`.

        Returns
        -------
        list
            A list of scaffold child nodes.

        """
        return self._get_scaffold_hierarchy(scaffold_smiles, data, default, max_levels, 'child')

    def add_scaffold_molecule_count(self):
        """Add the number of molecules containing each scaffold node as a scaffold
         node attribute ('count')."""
        for scaffold, data in self.nodes(data=True):
            mols = self.get_molecules_for_scaffold(scaffold)
            data['count'] = len(mols)

    def separate_disconnected_components(self, sort=False):
        """Separate disconnected components into distinct ScaffoldGraph objects.

        Parameters
        ----------
        sort : bool, optional
            If True sort components in descending order according to the
            number of nodes in the subgraph. The default is False.

        Returns
        -------
        list
            A list of ScaffoldGraph objects.

        """
        components = []
        for c in nx.weakly_connected_components(self):
            components.append(self.subgraph(c).copy())
        if sort:
            return sorted(components, key=len, reverse=True)
        return components

    def add_molecule_node(self, molecule, **attr):
        """Add a molecule node to the graph.

        Parameters
        ----------
        molecule : rdkit.Chem.rdchem.Mol
            Molecule to add to the graph. The molecule
            must have the '_Name' property available.
        **attr : keyword arguments, optional
            Attributes to add to the node.

        Notes
        -----
        The SMILES string of the molecule is added as
        a node attribute with the key 'smiles'.

        """
        name = molecule.GetProp('_Name')
        default_attr = dict(type='molecule', smiles=MolToSmiles(molecule))
        default_attr.update(molecule.GetPropsAsDict())
        default_attr.update(attr)
        self.add_node(name, **default_attr)

    def add_scaffold_node(self, scaffold, **attr):
        """Add a scaffold node to the graph.

        Parameters
        ----------
        scaffold : scaffoldgraph.core.Scaffold
            Scaffold to add to the graph.
        **attr : keyword arguments, optional
            Attributes to add to the node.

        Notes
        -----
        The hierarchy of the scaffold is added as
        a node attribute with the key 'hierarchy'.

        """
        default_attr = dict(type='scaffold', hierarchy=scaffold.rings.count)
        default_attr.update(attr)
        self.add_node(scaffold.get_canonical_identifier(), **default_attr)

    def add_molecule_edge(self, molecule, scaffold, **attr):
        """Add a scaffold -> molecule edge.

        Parameters
        ----------
        molecule : rdkit.Chem.rdchem.Mol
            Molecule to connect with scaffold.
            The molecule must have the '_Name' property
            available.
        scaffold : scaffoldgraph.core.Scaffold
            Scaffold to connect with molecule
        **attr : keyword arguments, optional
            Attributes to add to the edge.

        """
        name = molecule.GetProp('_Name')
        default_attr = dict(type=0)
        default_attr.update(attr)
        self.add_edge(scaffold.get_canonical_identifier(), name, **default_attr)

    def add_scaffold_edge(self, parent, child, **attr):
        """Add a scaffold (parent) -> scaffold (child) edge.

        Parameters
        ----------
        parent : scaffoldgraph.core.Scaffold
            Molecule to connect with scaffold
        child : scaffoldgraph.core.Scaffold
            Scaffold to connect with molecule
        **attr : keyword arguments, optional
            Attributes to add to the edge.

        """
        default_attr = dict(type=1)
        default_attr.update(attr)
        self.add_edge(
            parent.get_canonical_identifier(),
            child.get_canonical_identifier(),
            **default_attr
        )

    @classmethod
    def from_sdf(cls, file_name, ring_cutoff=10, progress=False, annotate=True,
                 zipped=False, flatten_isotopes=False, keep_largest_fragment=False,
                 discharge_and_deradicalize=False, **kwargs):
        """Construct a ScaffoldGraph from an SDF file.

        Parameters
        ----------
        file_name : str
            File path to an SDF input.
        ring_cutoff : int, optional
            Ignore molecules with more rings than this cutoff to avoid extended
            calculation time. The default is 10.
        progress : bool, optional
            If True display a progress bar to monitor construction progress.
            The default is False.
        annotate : bool, optional
            If True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold). The default is True.
        flatten_isotopes : bool, optional
            If True remove specific isotopes when initializing the scaffold.
            The default is False.
        keep_largest_fragment : bool, optional
            If True when encountering molecules containing disconnected fragments
            initialize the scaffold from only the largest disconnected fragment.
            The default is False.
        discharge_and_deradicalize : bool, optional
            If True remove charges and radicals when initializing the scaffold.
            The default is False.
        zipped : bool, optional
            If True input file is compressed with gzip. The default is False.
        **kwargs : keyword arguments, optional
            Arguments to pass to the ScaffoldGraph initializer.

        """
        if zipped:
            sdf = gzip.open(file_name, 'rb')
        else:
            sdf = open(file_name, 'rb')
        supplier = read_sdf(sdf, requires_length=progress is True)
        instance = cls(**kwargs)
        init_args = dict(
            flatten_isotopes=flatten_isotopes,
            keep_largest=keep_largest_fragment,
            discharge=discharge_and_deradicalize,
            annotate=annotate,
        )
        instance._construct(supplier, init_args, ring_cutoff=ring_cutoff, progress=progress)
        sdf.close()
        return instance

    @classmethod
    def from_smiles_file(cls, file_name, delimiter=' ', smiles_column=0, name_column=1, header=False,
                         ring_cutoff=10, progress=False, annotate=True, flatten_isotopes=False,
                         keep_largest_fragment=False, discharge_and_deradicalize=False, **kwargs):
        """Construct a ScaffoldGraph from a SMILES file.

        Parameters
        ----------
        file_name : str
            File path to a SMILES file input.
        delimiter : str, optional
            Delimiter used in SMILES file. The default is ' '.
        smiles_column : int, optional
            Index of column containing SMILES strings. The default is 0.
        name_column : int, optional
            Index of column containing molecule names. The default is 1.
        header : bool, optional
            If True skip the first line of the SMILES file. The default is False.
        ring_cutoff : int, optional
            Ignore molecules with more rings than this cutoff. The default is 10.
        progress : bool, optional
            If True display a progress bar to monitor construction progress.
            The default is False.
        annotate : bool, optional
            If True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold). The default is True.
        flatten_isotopes : bool, optional
            If True remove specific isotopes when initializing the scaffold.
            The default is False.
        keep_largest_fragment : bool, optional
            If True when encountering molecules containing disconnected fragments
            initialize the scaffold from only the largest disconnected fragment.
            The default is False.
        discharge_and_deradicalize : bool, optional
            If True remove charges and radicals when initializing the scaffold.
            The default is False.
        **kwargs : keyword arguments, optional
            Arguments to pass to the ScaffoldGraph initializer.

        """
        supplier = read_smiles_file(
            file_name, delimiter, smiles_column, name_column,
            header, requires_length=progress is True
        )
        instance = cls(**kwargs)
        init_args = dict(
            flatten_isotopes=flatten_isotopes,
            keep_largest=keep_largest_fragment,
            discharge=discharge_and_deradicalize,
            annotate=annotate,
        )
        instance._construct(supplier, init_args, ring_cutoff=ring_cutoff, progress=progress)
        return instance

    @classmethod
    def from_supplier(cls, supplier, ring_cutoff=10, progress=False, annotate=True,
                      flatten_isotopes=False, keep_largest_fragment=False,
                      discharge_and_deradicalize=False, **kwargs):
        """Construct a ScaffoldGraph from a custom rdkit Mol supplier.

        A simple supplier could be a list of rdkit molecules or a supplier provided by rdkit.

        Parameters
        ----------
        supplier : iterable
            A custom supplier of rdkit molecules (rdkit.Chem.rdchem.Mol), it is expected that
            each molecule will be assigned a property '_Name' serving as an identifier for that
            molecule. To use the progress option the supplier must have an associated length.
        ring_cutoff : int, optional
            Ignore molecules with more rings than this cutoff. The default is 10.
        progress : bool, optional
            If True display a progress bar to monitor construction progress.
            The default is False.
        annotate : bool, optional
            If True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold). The default is True.
        flatten_isotopes : bool, optional
            If True remove specific isotopes when initializing the scaffold.
            The default is False.
        keep_largest_fragment : bool, optional
            If True when encountering molecules containing disconnected fragments
            initialize the scaffold from only the largest disconnected fragment.
            The default is False.
        discharge_and_deradicalize : bool, optional
            If True remove charges and radicals when initializing the scaffold.
            The default is False.
        **kwargs : keyword arguments, optional
            Arguments to pass to the ScaffoldGraph initializer.

        Notes
        -----
        If the supplied molecules do not contain the '_Name' attribute a hash string will
        be used as the molecule node key.

        """
        instance = cls(**kwargs)
        init_args = dict(
            flatten_isotopes=flatten_isotopes,
            keep_largest=keep_largest_fragment,
            discharge=discharge_and_deradicalize,
            annotate=annotate,
        )
        instance._construct(supplier, init_args, ring_cutoff=ring_cutoff, progress=progress)
        return instance

    @classmethod
    def from_dataframe(cls, df, smiles_column='Smiles', name_column='Name', data_columns=None,
                       ring_cutoff=10, progress=False, annotate=True, flatten_isotopes=False,
                       keep_largest_fragment=False, discharge_and_deradicalize=False, **kwargs):
        """Construct a ScaffoldGraph from a pandas DataFrame.

        Parameters
        ----------
        df : pd.DataFrame
            A pandas DataFrame containing SMILES strings and a molecule identifier.
        smiles_column : value, optional
            Label of column containing SMILES strings. The default is 'Smiles'.
        name_column : str
            Label of column containing SMILES strings. The default is 'Name'.
        data_columns : list
            List of column keys to be included in the molecule node attributes.
        ring_cutoff : int, optional
            Ignore molecules with more rings than this cutoff. The default is 10.
        progress : bool, optional
            If True display a progress bar to monitor construction progress.
            The default is False.
        annotate : bool, optional
            If True write an annotated murcko scaffold SMILES string to each
            molecule edge (molecule --> scaffold). The default is True.
        flatten_isotopes : bool, optional
            If True remove specific isotopes when initializing the scaffold.
            The default is False.
        keep_largest_fragment : bool, optional
            If True when encountering molecules containing disconnected fragments
            initialize the scaffold from only the largest disconnected fragment.
            The default is False.
        discharge_and_deradicalize : bool, optional
            If True remove charges and radicals when initializing the scaffold.
            The default is False.
        **kwargs : keyword arguments, optional
            Arguments to pass to the ScaffoldGraph initializer.

        """
        supplier = read_dataframe(df, smiles_column, name_column, data_columns)
        instance = cls(**kwargs)
        init_args = dict(
            flatten_isotopes=flatten_isotopes,
            keep_largest=keep_largest_fragment,
            discharge=discharge_and_deradicalize,
            annotate=annotate,
        )
        instance._construct(supplier, init_args, ring_cutoff=ring_cutoff, progress=progress)
        return instance

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )
