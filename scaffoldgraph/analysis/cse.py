"""
scaffoldgraph.analysis.cse
"""

from concurrent.futures import ProcessPoolExecutor
from functools import partial
from itertools import repeat

import networkx as nx
import pandas as pd
import tqdm
from scipy.stats import ks_2samp, binom_test

from ..core.graph import ScaffoldGraph


class CSE(object):
    """Compound set enrichment (CSE)"""

    def __init__(self, hypothesis_test='ks', bonferroni=True,
                 p=0.5, crit=0.05, alternative='greater', n_jobs=1):

        """ Initialize CSE.

        Parameters
        ----------
        hypothesis_test : (str)
            if 'ks' use Kolmogorov-Smirnoff if 'binom' use a binomial test
        bonferroni : (bool)
            if True update the critical value for each scaffold using bonferroni
            correction applied to each hierarchy of the graph
        p : (float)
            The hypothesized probability of success. 0 <= p <= 1 for binomial test
            default: 0.5
        crit : (float)
            The critical level of significance default: 0.05
        alternative : (str) {'greater', 'less', 'two-sided'}
            Indicates the alternative hypothesis. default: 'greater'
        n_jobs : (int)
            Number of processes to use in calculation. default: 1
        """

        self._hyp_str = hypothesis_test
        if hypothesis_test == 'ks':
            self._hyp = ks_2samp
        elif hypothesis_test == 'binom':
            self._hyp = binom_test
        else:
            raise ValueError(f'hypothesis test {hypothesis_test}'
                             f' not currently implemented')

        alternatives = ['greater', 'less', 'two-sided']
        if alternative not in alternatives:
            raise ValueError(f'alternative must be one of: {alternatives}')
        self._alt = alternative

        self._bonferroni = bonferroni
        self._crit = crit
        self._p = p

        self._n_jobs = n_jobs

        self.is_fit = False
        self._X = dict()

    def fit(self, graph, bioactivity_map, progress=True):
        """Fit CSE to a scaffold graph.

        Parameters
        ----------
        graph : (sg.ScaffoldGraph)
            constructed scaffold graph
        bioactivity_map : (dict {str: float})
            dictionary containing a mapping of molecule IDs to
            a bio-activity value i.e. {MOLID: bioactivity}
        progress : (bool)
            if True show a progress bar to monitor progress
        """

        if not issubclass(type(graph), ScaffoldGraph):
            raise ValueError('Input graph must be a subclass of ScaffoldGraph')

        mapping = {k: d for k, d in bioactivity_map.items() if k in graph.nodes}

        def get_activity(scaffold):
            distribution = []
            succ = nx.bfs_tree(graph, scaffold)
            for node in succ.nodes:
                if graph.nodes[node]['type'] == 'molecule':
                    try:
                        distribution.append(mapping[node])
                    except KeyError:
                        continue
            return scaffold, distribution

        activity = (get_activity(x) for x in graph.get_scaffold_nodes())

        if self._hyp_str == 'ks':
            self._fit_ks(list(mapping.values()), activity, progress)
        elif self._hyp_str == 'binom':
            self._fit_binom(activity, progress)

        if self._bonferroni:
            hier = graph.get_hierarchy_sizes()
            for x in hier:
                hier[x] = self._crit / hier[x]
            for x in self._X.keys():
                h = graph.nodes[x]['hierarchy']
                self._X[x]['CRIT'] = hier[h]
        else:
            for x in self._X.keys():
                self._X[x]['CRIT'] = self._crit

        self.is_fit = True

    def _fit_ks(self, background, activity, p):
        scaffolds, activity = zip(*activity)
        func = partial(self._hyp, alternative=self._alt)
        with ProcessPoolExecutor(self._n_jobs) as pool:
            for f in tqdm.tqdm(zip(pool.map(func, repeat(background), activity), scaffolds),
                               total=len(scaffolds), disable=p is False):
                k, p = f[0]
                self._X.update({f[1]: {'PVAL': p, 'KS': k}})

    def _fit_binom(self, activity, p):
        activity = ((x[0], [x[1].count(1), x[1].count(0)]) for x in activity)
        scaffolds, activity = zip(*activity)
        func = partial(self._hyp, n=None, p=self._p, alternative=self._alt)
        with ProcessPoolExecutor(self._n_jobs) as pool:
            for f in tqdm.tqdm(zip(pool.map(func, activity), scaffolds),
                               total=len(scaffolds), disable=p is False):
                self._X.update({f[1]: {'PVAL': f[0]}})

    def to_dataframe(self, top_n=100):
        """Write CSE results to a dataframe.

        dataframe contains three columns:
            SCAFFOLD: scaffold canonical identifier
            PVAL: p-value
            CRIT: critical level of significance

        Parameters
        ----------
        top_n : (int)
            write the top N results to the dataframe starting
            with the lowest p-value.
        """

        if not self.is_fit:
            raise ValueError('CSE is not fit')

        df = pd.DataFrame()
        x = sorted(self._X, key=lambda x: self._X[x]['PVAL'])
        scaffolds = x[:top_n] if len(x) >= top_n else x
        df['SCAFFOLD'] = scaffolds
        df['PVAL'] = [self._X[s]['PVAL'] for s in scaffolds]
        df['CRIT'] = [self._X[s]['CRIT'] for s in scaffolds]

        return df

    def __getitem__(self, item):
        return self._X[item]

    def __len__(self):
        return len(self._X.keys())

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )
