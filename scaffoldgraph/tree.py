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

    References
    ----------
    .. [1] Schuffenhauer, A., Ertl, P., Roggo, S., Wetzel, S., Koch, M. A., and Waldmann, H. (2007).
           The scaffold tree visualization of the scaffold universe by hierarchical scaffold classification.
           Journal of Chemical Information and Modeling, 47(1), 47–58. PMID: 17238248.
    """

    def __init__(self, graph=None, prioritization_rules=None):
        super(ScaffoldTree, self).__init__(graph, MurckoRingFragmenter(True))
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
        """Return the prioritization ruleset used"""
        return self.rules


def tree_frags_from_mol(mol, prioritization_rules=None):
    """Generate a scaffold tree from a single molecule without using networkx

    Parameters
    ----------
    mol: rdkit molecule for processing
    prioritization_rules: rules for prioritizing parent scaffolds
        (optional, default: None)

    Returns
    -------
    parents: an ordered list of rdkit Mols representing a scaffold tree
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
