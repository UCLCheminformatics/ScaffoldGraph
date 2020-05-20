"""
scaffoldgraph.analysis.representation

Module contains an adaption of the Automated Identification of Over-Represented
Scaffold Classes in HTS Data method from the paper: 'HierS: Hierarchical Scaffold Clustering
Using Topological Chemical Graphs'
"""

from networkx import set_node_attributes
from collections import OrderedDict
from itertools import combinations

from rdkit import DataStructs
from rdkit import Chem


class Cache(OrderedDict):
    """A basic implementation of an LRU cache using OrderedDict"""

    def __init__(self, maxsize=None, *args, **kwargs):
        self._maxsize = maxsize
        super(Cache, self).__init__(*args, **kwargs)

    @property
    def maxsize(self):
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

    def __repr__(self):
        return '{}(maxsize={})'.format(
            self.__class__.__name__,
            self.maxsize
        )


class MolecularSimilarityCache(object):

    __slots__ = ('_fp_func', '_sim_func', '_fp_cache', '_sim_cache')

    def __init__(self, fp_func=None, sim_func=None, fp_cache_maxsize=None, sim_cache_maxsize=None):
        self._fp_func = fp_func if fp_func else Chem.RDKFingerprint
        assert callable(self._fp_func), 'fp_func must be callable or None'
        self._sim_func = sim_func if sim_func else DataStructs.TanimotoSimilarity
        assert callable(self._sim_func), 'sim_func must be callable or None'
        self._fp_cache = Cache(fp_cache_maxsize)
        self._sim_cache = Cache(sim_cache_maxsize)

    @property
    def fp_func(self):
        return self._fp_func

    @fp_func.setter
    def fp_func(self, fp_func):
        setattr(self, '_fp_func', fp_func)
        self.clear()  # clear both caches

    @property
    def sim_func(self):
        return self._sim_func

    @sim_func.setter
    def sim_func(self, sim_func):
        setattr(self, '_sim_func', sim_func)
        self.clear_sim_cache()  # clear only similarity cache

    def get_fingerprint(self, mol_node):
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
        id1, id2 = mol_node_1[0], mol_node_2[0]
        key = tuple(sorted([id1, id2]))
        if key in self._sim_cache:
            return self._sim_cache[key]
        fp1 = self.get_fingerprint(mol_node_1)
        fp2 = self.get_fingerprint(mol_node_2)
        sim = self._sim_cache.setdefault(key, self.sim_func(fp1, fp2))
        return sim

    def clear_fp_cache(self):
        self._fp_cache.clear()

    def clear_sim_cache(self):
        self._sim_cache.clear()

    def clear(self):
        self.clear_fp_cache()
        self.clear_sim_cache()

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
    In this implementation it is knows as APS as the function enables the user to specify
    similarity metrics other than Tanimoto using the sim_func argument.

    Parameters
    ----------
    scaffoldgraph : (ScaffoldGraph)
    fp_func : (callable or None, optional (default=None)) A callable returning a
        molecular fingerprint from an RDKit Mol object. If None the fingerprint
        is an RDKFingerprint with default parameters.
    sim_func : (callable or None, optional (default=None)) A callable returning a
        similarity value (float) for a pair of fingerprint objects calculated by
        fp_func. If None the default metric is Tanimoto.
    skip_levels : (tuple, list, optional (default=None)) Skip any scaffolds in hierarchy
        levels specified. The aps and membership is set to 0.
    fp_cache_maxsize : (int or None, optional (default=None)) If a value is
        specified the maximum number of fingerprints cached at any time is equal
        to value. If set to None no limit is set.
    sim_cache_maxsize : (int or None, optional (default=None)) If a value is
        specified the maximum number of pairwise similarities cached at any time is
        equal to value. If set to None no limit is set.

    Returns
    -------
    A dict of dicts in the format {scaffold: {members, aps}} where members is the
    number of molecules in the scaffold cluster and aps is the average pairwise
    similarity of the molecules in the cluster.
    """
    aps_dict = {}
    cache = MolecularSimilarityCache(fp_func, sim_func, fp_cache_maxsize, sim_cache_maxsize)

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

    cache.clear()  # empty the cache
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
    scaffoldgraph : (ScaffoldGraph)
    threshold : (float, optional (default=0.80))
        Similarity threshold used to define potential
        over-represented scaffolds.
    member_cutoff : (int or None, optional (default=None))
        If set, scaffolds for which (member_cutoff <= member molecules) are not considered
        to be over-represented (not significant).
    skip_aps : (bool, optional (default=False))
        If True function assumes that the APS has already been calculated and 'members' and
        'aps' are scaffold node attributes (i.e. use if running the same function more than once
        with different thresholds).
    **kwargs :
        Arguments for the calc_average_pairwise_similarity function (calculating the APS metric)

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
