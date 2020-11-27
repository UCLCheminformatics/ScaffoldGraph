"""
scaffoldgraph.vis.base
"""

import networkx as nx

from abc import ABC, abstractmethod

from scaffoldgraph.core import ScaffoldGraph
from scaffoldgraph.utils import canonize_smiles


class Visualizer(ABC):
    """Base class for ScaffoldGraph visualizers.

    A Visualizer contains functions for creating visualizations
    of ScaffoldGraphs.

    See Also
    --------
    scaffoldgraph.vis.notebook.cytoscape.CytoscapeVisualizer

    """
    def __init__(self, graph, requires_tree=False):
        """Initialize the visualizer.

        Parameters
        ----------
        graph : ScaffoldGraph
            ScaffoldGraph to visualize
        requires_tree : bool, optional
            Whether the visualizer requires a tree
            structure to create a visualization.

        """
        self._requires_tree = requires_tree
        self._graph = self._validate_graph(graph)

    @property
    def graph(self):
        """ScaffoldGraph: return the graph associated with the visualizer."""
        return self.graph

    @graph.setter
    def graph(self, graph):
        self._graph = self._validate_graph(graph)

    def _validate_graph(self, graph):
        """Private: Validate a graph is suitable for visualizer."""
        if not issubclass(type(graph), ScaffoldGraph):
            raise ValueError(
                f'{graph} must be a subclass of ScaffoldGraph'
            )
        if self._requires_tree:
            if not nx.is_tree(graph) or nx.is_forest(graph):
                msg = '{} requires a tree/forest structured graph'
                msg.format(self.__class__.__name__)
                raise ValueError(msg)
        return graph

    def _subgraph_from_mol(self, molecule):
        """Private: Select a subgraph starting at a molecule node.

        Parameters
        ----------
        molecule : str
            Molecule node identifier.

        Returns
        -------
        subgraph : ScaffoldGraph
            A subgraph starting at `molecule`.

        """
        G = self._graph
        if not G.molecule_in_graph(molecule):
            raise ValueError(f'molecule: {molecule} not in graph {G}')
        scaffolds = G.get_scaffolds_for_molecule(molecule)
        subgraph = G.subgraph([molecule] + scaffolds)
        return subgraph

    def _subgraph_from_scf(self, scaffold, traversal):
        """Private: Select a subgraph starting at a scaffold node.

        Parameters
        ----------
        scaffold : str
            Scaffold node identifier.
        traversal : str {'parent', 'child', 'bidirectional'}
            The direction of traversal to create the subgraph.
            If 'bidirectional' both directions are considered.

        Returns
        -------
        subgraph : ScaffoldGraph
            A subgraph starting at `scaffold`.

        """
        G = self._graph
        query = canonize_smiles(scaffold)
        if not G.scaffold_in_graph(query):
            raise ValueError(f'scaffold: {query} not in graph {G}')
        if traversal == 'parent':
            nodes = G.get_parent_scaffolds(query)
        elif traversal == 'child':
            nodes = list(nx.descendants(G, query))
        elif traversal == 'bidirectional':
            nodes = G.get_parent_scaffolds(query)
            nodes += list(nx.descendants(G, query))
        else:
            msg = 'traversal must be one of {child, parent, bidirectional}'
            raise ValueError(msg)
        subgraph = G.subgraph([query] + nodes)
        return subgraph

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )
