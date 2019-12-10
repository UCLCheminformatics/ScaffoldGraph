"""
scaffoldgraph.network
"""

from .core import MurckoRingFragmenter, MurckoRingSystemFragmenter
from .core import ScaffoldGraph


class ScaffoldNetwork(ScaffoldGraph):
    """
    Class representing a scaffold network.

    References
    ----------
    .. [1] Varin, T., Schuffenhauer, A., Ertl, P., and Renner, S. (2011). Mining for bioactive
           scaffolds with scaffold networks: Improved compound set enrichment from primary screening data.
           Journal of Chemical Information and Modeling, 51(7), 1528–1538.
    """

    def __init__(self, graph=None):
        super(ScaffoldNetwork, self).__init__(graph, MurckoRingFragmenter())

    def _recursive_constructor(self, child):
        parents = (p for p in self.fragmenter.fragment(child) if p)
        for parent in parents:
            if parent in self.nodes:
                self.add_scaffold_edge(parent, child)
            else:
                self.add_scaffold_node(parent)
                self.add_scaffold_edge(parent, child)
                if parent.rings.count > 1:
                    self._recursive_constructor(parent)


class HierS(ScaffoldGraph):
    """
    Class representing a HierS type scaffold network.

    Notes
    -----
    A HierS type network differs from a conventional scaffold network, through construction.
    When fragmenting molecules the HierS constructor does not attempt to break fused ring
    systems.

    References
    ----------
    .. [1] Wilkens, S., Janes, J., and Su, A. (2005). HierS: Hierarchical Scaffold Clustering
           Using Topological Chemical Graphs. Journal of Medicinal Chemistry, 48(9), 3182-3193.
    """

    def __init__(self, graph=None):
        super(HierS, self).__init__(graph, MurckoRingSystemFragmenter())

    def _recursive_constructor(self, child):
        parents = (p for p in self.fragmenter.fragment(child) if p)
        for parent in parents:
            if parent in self.nodes:
                self.add_scaffold_edge(parent, child)
            else:
                self.add_scaffold_node(parent)
                self.add_scaffold_edge(parent, child)
                if parent.ring_systems.count > 1:
                    self._recursive_constructor(parent)
