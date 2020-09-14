"""
scaffoldgraph.analysis.representation

Module contains an adaption of the Automated Identification of Over-Represented
Scaffold Classes in HTS Data method from the paper: 'HierS: Hierarchical Scaffold Clustering
Using Topological Chemical Graphs'
"""

from networkx import set_node_attributes
from collections import OrderedDict
from itertools import combinations
from operator import eq as _eq

from rdkit import DataStructs
from rdkit import Chem


class Cache(OrderedDict):
    """A basic implementation of an LRU cache using OrderedDict.

    Adapted (slightly) from the collections ``OrderedDict``
    documentation.

    .. _collections OrderedDict Documentation:
   https://docs.python.org/3/library/collections.html#collections.OrderedDict

    """
    def __init__(self, maxsize=None, *args, **kwargs):
        """
        Parameters
        ----------
        maxsize : int, None, optional
            Set the maximum size of the cache, if None the cache
            has no size limitation. The default is None.
        *args
            Variable length argument list.
            Passed to OrderedDict.
        **kwargs
            Arbitrary keyword arguments.
            Passed to OrderedDict.

        """
        self._maxsize = maxsize
        super(Cache, self).__init__(*args, **kwargs)

    @property
    def maxsize(self):
        """int: The maximum size of the cache."""
        return self._maxsize

    def __getitem__(self, key):
        value = super().__getitem__(key)
        self.move_to_end(key)
        return value

    def __setitem__(self, key, value):
        super().__setitem__(key, value)
        if self.maxsize and len(self) > self.maxsize:
            oldest = next(iter(self))
            del self[oldest]

    def __eq__(self, other):
        if isinstance(other, Cache):
            return dict.__eq__(self, other) and all(map(_eq, self, other))
        return dict.__eq__(self, other)

    def __repr__(self):
        return '{}(maxsize={})'.format(
            self.__class__.__name__,
            self.maxsize
        )


class MolecularSimilarityCache(object):
    """An LRU cache for speeding up repeated molecular similarity computations."""

    __slots__ = ('_fp_func', '_sim_func', '_fp_cache', '_sim_cache')

    def __init__(self, fp_func=None, sim_func=None, fp_cache_maxsize=None, sim_cache_maxsize=None):
        """
        Parameters
        ----------
        fp_func : callable, optional
            A function calculating a molecular fingerprint from an rdkit Mol object.
            If None the function is set to ``rdkit.Chem.RDKFingerprint``.
        sim_func : callable, optional
            A function calculating the similarity between two fingerprints as returned
            by `fp_func`. If None the function is set to ``rdkit.Datastructs.TanimotoSimilarity``
        fp_cache_maxsize : int, optional
            Set the maximum number of fingerprints cached. If None the cache is unbounded.
        sim_cache_maxsize : int, optional
             Set the maximum number of similarity values cached. If None the cache is unbounded.

        """
        self._fp_func = fp_func if fp_func else Chem.RDKFingerprint
        assert callable(self._fp_func), 'fp_func must be callable or None'
        self._sim_func = sim_func if sim_func else DataStructs.TanimotoSimilarity
        assert callable(self._sim_func), 'sim_func must be callable or None'
        self._fp_cache = Cache(fp_cache_maxsize)
        self._sim_cache = Cache(sim_cache_maxsize)

    @property
    def fp_func(self):
        """callable: The fingerprinting function

        If the fingerprinting function is changed both the similarity and
        fingerpint caches are cleared.
        """
        return self._fp_func

    @fp_func.setter
    def fp_func(self, fp_func):
        setattr(self, '_fp_func', fp_func)
        self.clear()  # clear both caches

    @property
    def sim_func(self):
        """callable: The molecular similarity function

        If the similarity function is changed the similarity cache is cleared.
        """
        return self._sim_func

    @sim_func.setter
    def sim_func(self, sim_func):
        setattr(self, '_sim_func', sim_func)
        self.clear_sim_cache()  # clear only similarity cache

    def get_fingerprint(self, mol_node):
        """Retrieve a fingerprint from the cache if it exists else calculate.

        Parameters
        ----------
        mol_node : tuple
            A molecule node from a ScaffoldGraph where the first entry is the
            molecule ID and the second is a dictionary of node attributes.

        Returns
        -------
        object
            A molecular fingerprint.

        """
        mol_id = mol_node[0]
        if mol_id in self._fp_cache:
            return self._fp_cache[mol_id]
        smi = mol_node[1]['smiles']
        fp = self._fp_cache.setdefault(mol_id, self._calc_fp(smi))
        return fp

    def _calc_fp(self, smiles):
        rdmol = Chem.MolFromSmiles(smiles)
        return self._fp_func(rdmol)

    def get_similarity(self, mol_node_1, mol_node_2):
        """Retrieve a similarity value from the cache if it exists else calculate.

        Parameters
        ----------
        mol_node_1 : tuple
            A molecule node from a ScaffoldGraph where the first entry is the
            molecule ID and the second is a dictionary of node attributes.
        mol_node_2 : tuple
            A molecule node from a ScaffoldGraph where the first entry is the
            molecule ID and the second is a dictionary of node attributes.

        Returns
        -------
        float
            A molecular similarity score.

        """
        id1, id2 = mol_node_1[0], mol_node_2[0]
        key = tuple(sorted([id1, id2]))
        if key in self._sim_cache:
            return self._sim_cache[key]
        fp1 = self.get_fingerprint(mol_node_1)
        fp2 = self.get_fingerprint(mol_node_2)
        sim = self._sim_cache.setdefault(key, self.sim_func(fp1, fp2))
        return sim

    def clear_fp_cache(self):
        """Empty the fingerprint cache."""
        self._fp_cache.clear()

    def clear_sim_cache(self):
        """Empty the similarity cache."""
        self._sim_cache.clear()

    def clear(self):
        """Empty both the fingerprint and similarity caches."""
        self.clear_fp_cache()
        self.clear_sim_cache()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.clear()

    def __repr__(self):
        return '{}({}, {})'.format(
            self.__class__.__name__,
            repr(self._fp_cache),
            repr(self._sim_cache)
        )


def calc_average_pairwise_similarity(scaffoldgraph, fp_func=None, sim_func=None, skip_levels=None,
                                     fp_cache_maxsize=None, sim_cache_maxsize=None):

    """Calculate average pairwise similarity for each scaffold in a ScaffoldGraph.

    Average Pairwise Similarity (APS) is a simple method for approximating the overall topological
    similarity between compounds in a given scaffold class. The APS coefficient can also be used
    as a metric to gauge scaffold over-representation in a set of compounds as described in the
    HierS paper.

    Notes
    -----
    The metric used in the HierS implementation is called APT (Average Pairwise Tanimoto).
    In this implementation it is known as 'APS', as the function enables the user to specify
    similarity metrics other than Tanimoto using the `sim_func` argument.

    Parameters
    ----------
    scaffoldgraph : ScaffoldGraph
    fp_func : callable, None, optional
        A callable returning a molecular fingerprint from an RDKit Mol object.
        If None the fingerprint is an RDKFingerprint with default parameters.
    sim_func : callable, None, optional
        A callable returning a similarity value (float) for a pair of fingerprint objects
        calculated by `fp_func`. If None the default metric is Tanimoto.
    skip_levels : iterable, None, optional
        Skip any scaffolds in hierarchy levels specified.
        The aps and membership is set to 0.
    fp_cache_maxsize : int, optional
            Set the maximum number of fingerprints cached. If None the cache is unbounded.
    sim_cache_maxsize : int, optional
             Set the maximum number of similarity values cached. If None the cache is unbounded.

    Returns
    -------
    dict
        A dict of dicts in the format {scaffold: {members, aps}} where members is the
        number of molecules in the scaffold cluster and aps is the average pairwise
        similarity of the molecules in the cluster.

    See Also
    --------
    scaffoldgraph.analysis.representation.get_over_represented_scaffold_classes

    """
    aps_dict = {}
    cache_args = (fp_func, sim_func, fp_cache_maxsize, sim_cache_maxsize)

    with MolecularSimilarityCache(*cache_args) as cache:
        for scaffold, data in scaffoldgraph.get_scaffold_nodes(True):
            aps_data = aps_dict.setdefault(scaffold, {})

            if skip_levels and data['hierarchy'] in skip_levels:
                aps_data['members'] = 0
                aps_data['aps'] = 0.0

            m_nodes = scaffoldgraph.get_molecules_for_scaffold(scaffold, data=True)
            n_members = len(m_nodes)
            aps_data['members'] = n_members

            # If only 1 member (or less in case of disconnect) set aps to 0.0
            if n_members <= 1:
                aps_data['aps'] = 0.0
                continue

            pw_sims = []
            for i, j in combinations(m_nodes, 2):
                pw_sims.append(cache.get_similarity(i, j))
            aps_data['aps'] = sum(pw_sims) / len(pw_sims)

    return aps_dict


def get_over_represented_scaffold_classes(scaffoldgraph, threshold=0.80, member_cutoff=None,
                                          skip_aps=False, **kwargs):

    """Returns scaffolds that are potentially over-represented in the dataset.

    This method is an adaptation of the method described in the HierS paper for
    automated identification of over-represented scaffold classes in HTS data.

    The algorithm first builds a list of all scaffolds exceeding the user-defined
    similarity threshold which is subsequently ordered by ascending scaffold hierarchy (HierS
    used molecular weight to sort, but using hierarchy makes sense as it is pre-calculated
    during construction). Each scaffold (above hierarchy 1) is then inspected to see if it is
    derived from any scaffold that precedes it in the list. Any scaffold in the list of
    overrepresented scaffolds that is found to be derived from a higher ranking (i.e.,
    lower molecular weight) scaffold is removed because all of the compounds that have membership
    in such scaffolds are already accounted for by the higher ranking scaffold.

    The HierS paper uses three defined similarity thresholds (APS) in three categories:

        loose  = 0.75
        medium = 0.80
        strict = 0.85

    Parameters
    ----------
    scaffoldgraph : ScaffoldGraph
    threshold : float, optional
        Similarity threshold used to define potential over-represented scaffolds.
        The default is 0.80 (i.e. medium)
    member_cutoff : int, None, optional
        If set, scaffolds for which (member_cutoff <= member molecules) are not considered
        to be over-represented (not significant). The default is None.
    skip_aps : bool, optional
        If True the function assumes that the APS has already been calculated and 'members' and
        'aps' are scaffold node attributes (i.e. use if running the same function more than
         once with different thresholds). The default is False.
    **kwargs :
        Arguments for the calc_average_pairwise_similarity function (calculating the APS metric).

    References
    ----------
    .. [1] Wilkens, S., Janes, J., and Su, A. (2005). HierS: Hierarchical Scaffold Clustering
           Using Topological Chemical Graphs. Journal of Medicinal Chemistry, 48(9), 3182-3193.

    """
    if skip_aps is False:
        aps = calc_average_pairwise_similarity(scaffoldgraph, **kwargs)
        set_node_attributes(scaffoldgraph, aps)
        aps.clear()

    or_scaffolds = []
    for scaffold, d in scaffoldgraph.get_scaffold_nodes(data=True):
        if d.get('aps', 0) > threshold and not (member_cutoff and not d.get('members') >= member_cutoff):
            or_scaffolds.append((scaffold, d))
    or_scaffolds.sort(key=lambda n: n[1].get('hierarchy'))
    or_set = set([s for s, _ in or_scaffolds])

    def _filter(scaffold):
        s, data = scaffold
        if data.get('hierarchy', 1) == 1:
            return True
        elif any([p in or_set for p in scaffoldgraph.get_parent_scaffolds(s)]):
            return False
        return True

    return tuple(filter(_filter, or_scaffolds))
