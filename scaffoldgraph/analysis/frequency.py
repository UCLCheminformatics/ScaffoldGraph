"""
scaffoldgraph.analysis.frequency

"""
import numpy as np


def cumulative_scaffold_frequency(
        scaffoldgraph,
        hierarchy=-1,
        norm_hierarchy=False,
        frequency_key=None
):
    """Calculate cumulative scaffold frequency distrubutions (CSF) from
    a scaffold graph.

    Parameters
    ----------
    scaffoldgraph : ScaffoldGraph
        A ScaffoldGraph object to query.
    hierarchy : int
        The scaffold hierarchy to consider. If -1 then the CSF is
        calculated for murcko scaffolds rather than a scaffold
        hierarchy. The default is -1 (murcko scaffolds).
    norm_hierarchy : bool
        Normalise the CSF by the number of molecules represented
        by the considered hierarchy rather than the total molecules
        in the graph. If False then compound representation in the CDF
        may be below 100%. The default is False (normalise by total).
    frequency_key : str, None, optional
        If scaffold frequency exists as an attribute of the graph, set
        this key to avoid re-calculation of scaffold frequencies.

    Examples
    --------
    Create a CSF plot (CSFP) for murcko scaffolds

    >>> import matplotlib.pyplot as plt
    >>> fig, ax = plt.subplots()
    >>> x, y = cumulative_scaffold_frequency(tree, hierarchy=-1)
    >>> ax.plot(x, y, label='Murcko CSF')
    >>> ax.set_xlabel('Percentage of scaffolds')
    >>> ax.set_ylabel('Percentage of molecules')
    >>> ax.legend()
    >>> fig.show()

    Calculate P_50, the percentage of scaffolds that represent 50% of compounds.

    >>> import numpy as np
    >>> p50 = np.interp(0.5, y, x)

    Returns
    -------
    tuple
        A tuple containing the cumulative percentage of scaffolds and
        the cumulative scaffold frequency as a percentage of molecules.
        Can be used to plot a CSF plot (x, y).

    Notes
    -----
    Cumulative scaffold frequency should be used with scaffold tree
    structures.

    """
    if not frequency_key:
        scaffoldgraph.add_scaffold_molecule_count()
        frequency_key = 'count'
    if hierarchy == -1:  # murcko scaffolds
        h = _get_murcko_frequency(scaffoldgraph)
    elif hierarchy in scaffoldgraph.get_hierarchy_sizes():  # hierarchies
        h = scaffoldgraph.get_scaffolds_in_hierarchy(hierarchy, frequency_key)
    else:  # hierarchy does not exist
        raise ValueError(f'Invalid hierarchy: {hierarchy}')
    sh = sorted(h, key=lambda x: x[1], reverse=True)
    sf = 1. * np.arange(len(sh)) / (len(sh) - 1)
    cumsum = np.cumsum([x[1] for x in sh])
    if norm_hierarchy is True and hierarchy > 0:
        mf = cumsum / cumsum[-1]
    else:  # normalize by total molecules in graph
        mf = cumsum / scaffoldgraph.num_molecule_nodes
    return sf, mf


def _get_murcko_frequency(scaffoldgraph):
    """Get frequencies for murcko scaffolds."""
    g = scaffoldgraph
    mols = g.get_molecule_nodes()
    m = list({next(g.predecessors(x)) for x in mols})
    f = [len([x for x in g.successors(x) if g.nodes[x]['type'] == 'molecule']) for x in m]
    return list(zip(m, f))
