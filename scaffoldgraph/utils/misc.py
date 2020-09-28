"""
scaffoldgraph.utils.misc

Defines miscellaneous functions used within scaffoldgraph.
"""

import networkx as nx

from rdkit import Chem


def canonize_smiles(smiles, failsafe=True):
    """Canonize a SMILES string (with failsafe).

    Parameters
    ----------
    smiles : str
        SMILES string to canonize.
    failsafe : bool
        If True, if the SMILES fails to parse
        the input SMILES is returned instead
        of raising an error.

    Returns
    -------
    str
        The canonical SMILES representation.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None and failsafe:
        return smiles
    return Chem.MolToSmiles(mol)


def summary(graph, n=None):
    """Return a summary of information for the graph or a single node n.

    Parameters
    ----------
    graph : sg.core.ScaffoldGraph or NetworkX graph
       A graph object which can either be a ScaffoldGraph graph or a NetworkX
       graph object.
    n : any hashable, optional
       A node in the graph. The default is None.

    Returns
    -------
    info : str
        A string containing the summary.

    Raises
    ------
    ValueError
        If n is not in the graph.

    """
    from scaffoldgraph.core import ScaffoldGraph
    if not issubclass(type(graph), ScaffoldGraph):
        return nx.info(graph, n)
    info = ""
    if n is None:
        type_name = [type(graph).__name__]
        info += f"Type: {','.join(type_name)}\n"
        info += f"Number of molecule nodes: {graph.num_molecule_nodes}\n"
        info += f"Number of scaffold nodes: {graph.num_scaffold_nodes}\n"
        info += f"Number of edges: {graph.number_of_edges()}\n"
        info += f"Max hierarchy: {graph.max_hierarchy()}\n"
        info += f"Min hierarchy: {graph.min_hierarchy()}\n"
    else:
        if graph.molecule_in_graph(n):
            info += f"Node {n} has the following properties:\n"
            info += "Type: molecule\n"
            info += f"SMILES: {graph.nodes[n].get('smiles')}\n"
            info += f"Degree: {graph.degree(n)}\n"
            info += "Parent scaffolds: "
            info += " ".join(str(s) for s in graph.predecessors(n))
        elif graph.scaffold_in_graph(n):
            key = canonize_smiles(n)
            info += f"Node {key} has the following properties:\n"
            info += "Type: scaffold\n"
            info += f"Hierarchy: {graph.nodes[key].get('hierarchy')}\n"
            info += f"Degree: {graph.degree(key)}\n"
            info += "Parent scaffolds: "
            info += " ".join(str(s) for s in graph.get_parent_scaffolds(key, max_levels=1))
            info += "\n"
            info += "Child scaffolds: "
            info += " ".join(str(s) for s in graph.get_child_scaffolds(key, max_levels=1))
            info += "\n"
            info += "Child molecules: "
            info += " ".join(
                str(s) for s in graph.successors(key) if graph.nodes[s].get('type') == 'molecule'
            )
        else:
            raise ValueError(f"node {n} not in graph")
    return info
