"""
scaffoldgraph.analysis.diversity
"""

from rdkit.SimDivFilters.rdSimDivPickers import LeaderPicker
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect
from rdkit.Chem import MolFromSmiles

from functools import partial


def _form_dist_func(dist_func, fps):
    """function: create a partial dist_func."""
    if dist_func.__code__.co_argcount != 3:
        raise ValueError('dist_func must have three arguments: i, j, fps')
    if dist_func.__code__.co_varnames[2] != 'fps':
        raise ValueError('dist_func third argument name must be: fps')
    formed_dist_func = partial(dist_func, fps=fps)
    return formed_dist_func


def _make_diversity_pick(pool, threshold, pick_size, dist_func=None):
    """iterable: make a diversity pick from a pool of fingerprints."""
    picker = LeaderPicker()
    pool_size = len(pool)
    if pick_size > pool_size:
        pick_size = pool_size
    if dist_func is None:
        pick = picker.LazyBitVectorPick(pool, pool_size, threshold, pick_size)
    else:
        dist_func = _form_dist_func(dist_func, pool)
        pick = picker.LazyPick(dist_func, pool_size, threshold, pick_size)
    return pick


def _create_pool(scaffold, graph, radius, bits):
    """tuple : create molecule pool (ids, mols, fps)."""
    mol_ids, smiles = zip(*graph.get_molecules_for_scaffold(scaffold, 'smiles'))
    mols = list(map(MolFromSmiles, smiles))
    fps = list(map(lambda x: GetMorganFingerprintAsBitVect(x, radius, nBits=bits), mols))
    if len(fps) == 0:
        raise ValueError(f'No molecules for scaffold class: {scaffold}')
    return mol_ids, mols, fps


def diversity_pick_for_scaffold_class(
        scaffold,
        graph,
        threshold=0.65,
        pick_size=0,
        fp_radius=2,
        fp_bits=1024,
        dist_func=None
):
    """
    Pick a diverse set of molecules from a scaffold class using
    the RDKit diversity picker (LeaderPicker) and Morgan
    fingerprints.

    Parameters
    ----------
    scaffold : str
        Scaffold class name i.e. scaffold SMILES.
    graph : ScaffoldGraph
        ScaffoldGraph for picking.
    threshold : float, optional
        Stop picking when the distance goes below this value.
        The default is 0.65 i.e. similarity = 0.35.
    pick_size : int, optional
        Number of items to pick from the molecule pool. If
        the pick size is greater than the pool size, the
        pick size will be equal to the size of the pool.
    fp_radius : int, optional
        Radius of Morgan fingerprint. The default is 2.
    fp_bits : int, optional
        Number of bits in the Morgan fingerprint. The
        default is 1024.
    dist_func : function, optional
        A function for calculating distance between a pair
        of fingerprints. The function should take two indicies
        (i, j) and a list of fingerprints (fps) and return
        the distance between these points.

    Examples
    --------
    Diversity pick for benzene scaffold.

    >>> ids, mols, fps = diversity_pick_for_scaffold_class('c1ccccc1', graph, pick_size=10)

    Returns
    -------
    tuple ((ids), (mols), (fps))
        A tuple of tuples with the first containg the picked molecules ids,
        the seconds containing the picked mols RDMols and the third containg
        the molecules fingerprints.

    Notes
    -----
    If performing diversity picks on a large scale, a custom implementation
    should probably be used where fingerprints can be cached.

    """
    mol_ids, mols, fps = _create_pool(scaffold, graph, fp_radius, fp_bits)
    pick = _make_diversity_pick(fps, threshold, pick_size, dist_func)
    picked = [(mol_ids[x], mols[x], fps[x]) for x in pick]
    ids, mols, fps = zip(*picked)
    return ids, mols, fps
