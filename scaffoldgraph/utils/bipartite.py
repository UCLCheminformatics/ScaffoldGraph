"""
scaffoldgraph.utils.bipartite

Defines functions for creating bipartite graphs from scaffold graphs.
"""

from scaffoldgraph.core import ScaffoldGraph


def make_bipartite_graph(graph):
    """Collapse a scaffold hierarchy into a bipartite representation.

    Scaffold --> Molecule

    The returned output will inherit the class of the input graph.

    Parameters
    ----------
    graph : sg.core.ScaffoldGraph
        A scaffold graph template for producing a bipaertite
        graph.

    Returns
    -------
    sg.core.ScaffoldGraph
        Bipartite scaffoldgraph where the scaffold hierarchy
        has been collapsed.

    """
    if not issubclass(type(graph), ScaffoldGraph):
        raise ValueError(f'{graph} must be a ScaffoldGraph')
    graph_type = type(graph)
    B = graph_type(None)
    for scf, sdata in graph.get_scaffold_nodes(True):
        B.add_node(scf, **sdata)
        for mol, mdata in graph.get_molecules_for_scaffold(scf, True):
            if not B.molecule_in_graph(mol):
                B.add_node(mol, **mdata)
            B.add_edge(scf, mol)
    return B
