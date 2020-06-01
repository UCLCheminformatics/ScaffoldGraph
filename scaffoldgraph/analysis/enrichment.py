"""
scaffoldgraph.analysis.enrichment

Module contains an implementation of Compound Set Enrichment from the papers:
- Compound Set Enrichment: A Novel Approach to Analysis of Primary HTS Data.
- Mining for bioactive scaffolds with scaffold networks: Improved compound set enrichment from primary screening data.
"""

from networkx import set_node_attributes
from scipy.stats import ks_2samp, binom_test
from loguru import logger


def _btp(scaffoldgraph, activity_key, alternative):
    result, active, total = {}, 0, 0
    for m, a in scaffoldgraph.get_molecule_nodes(activity_key):
        if a == 1:
            active += 1
        total += 1
    pd = active / total
    logger.debug(f'(BTP) Total: {total}, Active: {active}, pd: {pd}')
    for scaffold in scaffoldgraph.get_scaffold_nodes():
        mols, acts = zip(*scaffoldgraph.get_molecules_for_scaffold(scaffold, activity_key))
        N, K = len(mols), acts.count(1)
        pval = binom_test(K, N, pd, alternative=alternative)
        logger.debug(f'(BTP) {scaffold}, {K}, {N}, {pval}')
        result[scaffold] = {'pval': pval, '_active': K, '_total': N}
    return result


def _ksp(scaffoldgraph, activity_key, alternative):
    result, background = {}, []
    for _, activity in scaffoldgraph.get_molecule_nodes(activity_key):
        background.append(activity)
    for scaffold in scaffoldgraph.get_scaffold_nodes():
        mols, acts = zip(*scaffoldgraph.get_molecules_for_scaffold(scaffold, activity_key))
        N = len(mols)
        dmax, pval = ks_2samp(acts, background, alternative, 'auto')
        logger.debug(f'(KSP) {scaffold}, {N}, {dmax}, {pval}')
        result[scaffold] = {'pval': pval, 'dmax': dmax, '_total': N}
    return result


def bonferroni_correction(scaffoldgraph, crit):
    """Returns bonferroni corrected significance level for each hierarchy"""
    hier = scaffoldgraph.get_hierarchy_sizes()
    return {k: crit / v for k, v in hier.items()}


def calc_scaffold_enrichment(scaffoldgraph, activity, mode='ks', alternative='greater'):
    """
    Calculate scaffold enrichment using the Kolmogorov-Smirnov or binomal test.

    Parameters
    ----------
    scaffoldgraph : (ScaffoldGraph)
    activity : (string)
        A node attribute key corresponding to an activity value. If the test is binomial
        this value should be binary (0 or 1)
    mode : (string, optional (default='ks'))
        A string specifying the test to use to determine scaffold enrichment. 'ks': Kolmogorov-
        Smirnov, 'binomal': binomial test
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available (default is 'greater'):
          * 'two-sided'
          * 'less': one-sided
          * 'greater': one-sided

    Returns
    -------
    A dict of dicts in the format {scaffold: {pval}} where pval is the calculated p-value
    """

    if mode == 'binomial' or mode == 'b':
        return _btp(scaffoldgraph, activity, alternative)
    elif mode == 'ks' or mode == 'k':
        return _ksp(scaffoldgraph, activity, alternative)
    else:
        raise ValueError(f'scaffold enrichment mode: {mode}, not implemented')


def compound_set_enrichment(scaffoldgraph, activity, mode='ks', alternative='greater', crit=0.01):
    """
    Perform compound set enrichment (CSE), calculating scaffolds enriched for bioactivity.

    Parameters
    ----------
    scaffoldgraph : (ScaffoldGraph)
    activity : (string)
        A node attribute key corresponding to an activity value. If the test is binomial
        this value should be binary (0 or 1)
    mode : (string, optional (default='ks'))
        A string specifying the test to use to determine scaffold enrichment. 'ks': Kolmogorov-
        Smirnov, 'binomal': binomial test
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available (default is 'greater'):
          * 'two-sided'
          * 'less': one-sided
          * 'greater': one-sided
    crit : (float, optional (default=0.01))
        The critical significance level

    Returns
    -------
    A tuple of 'enriched' scaffold classes

    References
    ----------
    .. [1] Varin, T., Schuffenhauer, A., Ertl, P., and Renner, S. (2011). Mining for bioactive scaffolds
           with scaffold networks: Improved compound set enrichment from primary screening data.
           Journal of Chemical Information and Modeling, 51(7), 1528–1538.
    .. [2] Varin, T., Gubler, H., Parker, C., Zhang, J., Raman, P., Ertl, P. and Schuffenhauer, A. (2010)
           Compound Set Enrichment: A Novel Approach to Analysis of Primary HTS Data.
           Journal of Chemical Information and Modeling, 50(12), 2067-2078.
    """
    set_node_attributes(scaffoldgraph, calc_scaffold_enrichment(scaffoldgraph, activity, mode, alternative))
    bonferroni = bonferroni_correction(scaffoldgraph, crit)
    result = []
    for scaffold, data in scaffoldgraph.get_scaffold_nodes(True):
        if data['pval'] < bonferroni[data['hierarchy']]:
            result.append((scaffold, data))
    return tuple(sorted(result, key=lambda x: x[1]['pval']))
