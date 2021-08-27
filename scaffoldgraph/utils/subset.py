"""
scaffoldgraph.utils.subset

"""
from networkx.algorithms.traversal import bfs_tree
from collections import defaultdict


def split_graph_by_molecule_attribute(graph, key, default=None):
    """Split a scaffold graph into subgraphs based on unique molecule attributes.

    This function first groups molecule nodes sharing a unique attribute
    value, and then proceeds to build subgraphs from each node subset using
    a breadth-first search.

    The returned subgraphs are graph views and thus changes to the graph are
    nruled out by the view, but changes to node attributes
    are reflected in the original graph. To prevent this behaviour use:
    subgraph.copy()

    Parameters
    ----------
    graph : sg.core.ScaffoldGraph
        A scaffold graph to split.
    key : str
        The key for the molecule node attribute used to split the graph
        into subgraphs.
    default : value, bool, optional
        Value used for nodes that don't have the requested attribute.

    Returns
    -------
    splits : dict
        A dictionary with keys representing unique node attributes and
        values representing the constructed subgraphs.

    """
    if isinstance(key, bool):
        raise ValueError('Attribute key cannot be a boolean type')
    splits = defaultdict(list)
    for node, attr in graph.get_molecule_nodes(key, default):
        splits[attr].append(node)
    splits.default_factory = None  # Not really required
    for attr, nodes in splits.items():
        bfs_subset = set()
        for node in nodes:
            bfs = bfs_tree(graph, node, reverse=True)
            bfs_subset.update(bfs)
        splits[attr] = graph.subgraph(bfs_subset)
    return splits
