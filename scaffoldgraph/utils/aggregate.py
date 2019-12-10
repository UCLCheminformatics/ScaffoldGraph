"""
scaffoldgraph.utils.aggregate
"""

import warnings

from networkx import Graph, compose


def aggregate(list_of_graphs):
    """Aggregate a list of graphs into one graph object.

    Graphs within the list must be a subclass of a networkx Graph object

    Parameters
    ----------
    list_of_graphs (list): A list of scaffold graphs for aggregation.

    Returns
    -------
    A graph type object with the same class as the first entry in
    the parameter list_of_graphs

    Raises
    ------
    ValueError:
        raises if an empty list is provided, instead of a list of graphs
    ValueError:
        raises if any entry in the list is not a subclass of nx.Graph

    Examples
    --------
    >>> g1 = sg.ScaffoldNetwork.from_sdf('g1.sdf')
    >>> print(g1.number_of_nodes())
    100
    >>> g2 = sg.ScaffoldNetwork.from_sdf('g2.sdf')
    >>> print(g2.number_of_nodes())
    50
    >>> g3 = sg.ScaffoldNetwork.from_sdf('g3.sdf')
    >>> print(g3.number_of_nodes())
    200
    >>> list_of_graphs = [g1, g2, g3]
    >>> aggregated_graph = aggregate(list_of_graphs)
    >>> print(aggregated_graph.number_of_nodes())
    325

    Notes
    -----
    The user is not prevented from aggregating multiple graphs of
    differing types, although this may lead to undesired behaviour.
    (i.e. aggregating a tree and a network is possible)

    Based on nx.compose_all:
    .. _Compose-all: https://networkx.github.io/documentation/stable/reference/algorithms/
    generated/networkx.algorithms.operators.all.compose_all.html
    """

    if not list_of_graphs:
        raise ValueError('Cannot apply aggregate to an empty list')
    graphs = iter(list_of_graphs)
    C = next(graphs)
    graph_type = type(C)
    for H in graphs:
        if not issubclass(type(H), Graph):
            raise ValueError('Can only aggregate graph type objects')
        if graph_type != type(H):
            warnings.warn('Attempting to aggregate graphs of different types '
                          f'({graph_type} & {type(H)}) '
                          'could result in undesired behaviour')
        C = compose(C, H)
    return C
