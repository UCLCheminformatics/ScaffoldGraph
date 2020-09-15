"""
scaffoldgraph.tree
"""

from rdkit import RDLogger
from rdkit.Chem import rdmolops

from .core import ScaffoldGraph, Scaffold, MurckoRingFragmenter
from .core.fragment import get_murcko_scaffold
from .prioritization import original_ruleset

rdlogger = RDLogger.logger()


class ScaffoldTree(ScaffoldGraph):
    """
    Class representing a scaffold tree.

    Explore scaffold-space through the iterative removal of the least-characteristic
    ring from a molecular scaffold. The output is a tree of molecular scaffolds.

    Examples
    --------
    Create a ScaffoldTree from a SMILES file.

    >>> import scaffoldgraph as sg
    >>> tree = sg.ScaffoldTree.from_smiles_file('my_file.smi', progress=True)
    >>> network.num_scaffold_nodes
    75

    Create a ScaffoldTree from an SDF.

    >>> tree = sg.ScaffoldTree.from_sdf('my_file.sdf', progress=True)

    If the SDF is zipped:

    >>> tree = sg.ScaffoldTree.from_sdf('my_file.sdf.gz', zipped=True)

    Get scaffold nodes:

    >>> list(tree.get_scaffold_nodes())
    ['O=C(OCOC(=O)c1cccc2ncn(Cc3ccccc3)c12)OC1CCCCC1',
    'O=C(OCOC(=O)c1cccc2nc[nH]c12)OC1CCCCC1',
    ...]

    Include node attributes:

    >>> list(tree.get_scaffold_nodes(data=True))
    [('O=C(OCOC(=O)c1cccc2ncn(Cc3ccccc3)c12)OC1CCCCC1', {'type': 'scaffold', 'hierarchy': 4}),
    ('O=C(OCOC(=O)c1cccc2nc[nH]c12)OC1CCCCC1', {'type': 'scaffold', 'hierarchy': 3}),
    ...]

    Get molecule nodes (use data=True to get attributes):

    >>> list(tree.get_molecule_nodes())
    ['DB00006',
     'DB00007',
     'DB00014',
     ...]

    References
    ----------
    .. [1] Schuffenhauer, A., Ertl, P., Roggo, S., Wetzel, S., Koch, M. A., and Waldmann, H. (2007).
           The scaffold tree visualization of the scaffold universe by hierarchical scaffold classification.
           Journal of Chemical Information and Modeling, 47(1), 47–58. PMID: 17238248.

    See Also
    --------
    ScaffoldGraph
    ScaffoldNetwork
    HierS

    """
    def __init__(self, graph=None, prioritization_rules=None, **kwargs):
        """Initialize a ScaffoldTree.

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
        prioritization_rules : ScaffoldRuleSet
            Ruleset for prioritizing parent scaffolds during tree
            construction.

        """
        super(ScaffoldTree, self).__init__(graph, MurckoRingFragmenter(True), 'tree')
        self.rules = prioritization_rules if prioritization_rules else original_ruleset

    def _recursive_constructor(self, child):
        parents = [p for p in self.fragmenter.fragment(child) if p]
        if not parents:
            return
        parent = self.rules(child, parents)
        if not parent:
            return
        deletion_rule = parent.prioritization_rule
        if parent in self.nodes:
            self.add_scaffold_edge(parent, child, rule=deletion_rule)
        else:
            self.add_scaffold_node(parent)
            self.add_scaffold_edge(parent, child, rule=deletion_rule)
            if parent.rings.count > 1:
                self._recursive_constructor(parent)

    @property
    def prioritization_rules(self):
        """ScaffoldRuleSet : Return the prioritization ruleset used."""
        return self.rules


def tree_frags_from_mol(mol, prioritization_rules=None):
    """Generate a scaffold tree from a single molecule without using networkx.

    Parameters
    ----------
    mol: rdkit.Chem.rdchem.Mol
        rdkit molecule for processing.
    prioritization_rules : ScaffoldRuleSet, optional
        rules for prioritizing parent scaffolds. If
        not supplied the original rules are used.
        The default is None.

    Returns
    -------
    parents
        An ordered list of rdkit Mols representing a scaffold tree.

    Examples
    --------
    Generating scaffold tree fragments:

    >>> from rdkit import Chem
    >>> smiles = 'Cc1[nH]cnc1Cn1cccc(-c2ccccc2O)c1=O'
    >>> molecule = Chem.MolFromSmiles(smiles)
    >>> frags = tree_frags_from_mol(molecule)

    """
    rdlogger.setLevel(4)
    scaffold = Scaffold(get_murcko_scaffold(mol))
    rdmolops.RemoveStereochemistry(scaffold.mol)
    parents = [scaffold]
    fragmenter = MurckoRingFragmenter(use_scheme_4=True)
    rules = prioritization_rules if prioritization_rules else original_ruleset

    def _next_scaffold(child):
        next_parents = [p for p in fragmenter.fragment(child) if p]
        if not next_parents:
            return
        next_parent = rules(child, next_parents)
        parents.append(next_parent)
        if next_parent.rings.count > 1:
            _next_scaffold(next_parent)

    _next_scaffold(scaffold)
    rdlogger.setLevel(3)

    return [p.mol for p in parents]
