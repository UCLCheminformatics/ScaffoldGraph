"""
scaffoldgraph.analysis.enrichment

Module contains an implementation of Compound Set Enrichment from the papers:
- Compound Set Enrichment: A Novel Approach to Analysis of Primary HTS Data.
- Mining for bioactive scaffolds with scaffold networks: Improved compound set enrichment from primary screening data.
"""

from networkx import set_node_attributes
from scipy.stats import ks_2samp, binom_test
from loguru import logger


def _btp(scaffoldgraph, activity_key, alternative, pd):
    """CSE - binomial test (used in cse functions)."""
    result, active, total = {}, 0, 0
    for m, a in scaffoldgraph.get_molecule_nodes(activity_key):
        if int(a) == 1:
            active += 1
        total += 1
    if pd is None:
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
    """CSE - Kolmogorov-Smirnov test (used in cse functions)."""
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
    """Returns bonferroni corrected significance level for each hierarchy.

    Parameters
    ----------
    scaffoldgraph : ScaffoldGraph
        A ScaffoldGraph object to query.
    crit : float
        The critical significance value to apply bonferroni correction at
        each scaffold hierarchy.

    Returns
    -------
    dict
        A dictionary containing the corrected critical significance value
        at each scaffold hierarchy {hierarchy: crit}.

    """
    hier = scaffoldgraph.get_hierarchy_sizes()
    return {k: crit / v for k, v in hier.items()}


def calc_scaffold_enrichment(scaffoldgraph, activity, mode='ks', alternative='greater', p=None):
    """
    Calculate scaffold enrichment using the Kolmogorov-Smirnov or binomal test.

    Parameters
    ----------
    scaffoldgraph : ScaffoldGraph
        A ScaffoldGraph object to query.
    activity : str
        A scaffold node attribute key corresponding to an activity value.
        If the test is binomial this value should be a binary attribute
        (0 or 1 / True or False).
    mode : {'ks', 'b'}, optional
        A string specifying the statistical test to perform. 'ks' specifies a
        Kolmogorov-Smirnov test and 'b' or 'binomial' specifies a binomial test.
        The default is 'ks'.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available:
          * 'two-sided'
          * 'less': one-sided
          * 'greater': one-sided
        The default is 'greater'.
    p : float, None, optional
        The hypothesized probability of success. 0 <= p <= 1. Used in binomial mode.
        If not specified p is set automatically (number of active / total compounds).
        The default is None.

    Returns
    -------
    dict
        A dict of dicts in the format {scaffold: {results}} where results is the set
        of results returned by the statistical test and scaffold is a scaffold node
        key corresponding to a scaffold in the ScaffoldGraph object.

    See Also
    --------
    scaffoldgraph.analysis.enrichment.compound_set_enrichment

    References
    ----------
    .. [1] Varin, T., Schuffenhauer, A., Ertl, P., and Renner, S. (2011). Mining for bioactive scaffolds
           with scaffold networks: Improved compound set enrichment from primary screening data.
           Journal of Chemical Information and Modeling, 51(7), 1528–1538.
    .. [2] Varin, T., Gubler, H., Parker, C., Zhang, J., Raman, P., Ertl, P. and Schuffenhauer, A. (2010)
           Compound Set Enrichment: A Novel Approach to Analysis of Primary HTS Data.
           Journal of Chemical Information and Modeling, 50(12), 2067-2078.

    """
    if mode == 'binomial' or mode == 'b':
        return _btp(scaffoldgraph, activity, alternative, p)
    elif mode == 'ks' or mode == 'k':
        return _ksp(scaffoldgraph, activity, alternative)
    else:
        raise ValueError(f'scaffold enrichment mode: {mode}, not implemented')


def compound_set_enrichment(scaffoldgraph, activity, mode='ks', alternative='greater', crit=0.01, p=None):
    """
    Perform compound set enrichment (CSE), calculating scaffolds enriched for bioactivity.

    Parameters
    ----------
    scaffoldgraph : ScaffoldGraph
        A ScaffoldGraph object to query.
    activity : str
        A scaffold node attribute key corresponding to an activity value.
        If the test is binomial this value should be a binary attribute
        (0 or 1 / True or False).
    mode : {'ks', 'b'}, optional
        A string specifying the statistical test to perform. 'ks' specifies a
        Kolmogorov-Smirnov test and 'b' or 'binomial' specifies a binomial test.
        The default is 'ks'.
    alternative : {'two-sided', 'less', 'greater'}, optional
        Defines the alternative hypothesis.
        The following options are available:
          * 'two-sided'
          * 'less': one-sided
          * 'greater': one-sided
        The default is 'greater'.
    crit : float, optional
        The critical significance level. The default is 0.01
    p : float, None, optional
        The hypothesized probability of success. 0 <= p <= 1. Used in binomial mode.
        If not specified p is set automatically (number of active / total compounds).
        The default is None.

    Returns
    -------
    A tuple of 'enriched' scaffold classes in the format: (scaffold, {data}) where data
    is the corresponding node attributes for the returned scaffold.

    Notes
    -----
    P-values are added as node attributes with the key 'pval'.

    References
    ----------
    .. [1] Varin, T., Schuffenhauer, A., Ertl, P., and Renner, S. (2011). Mining for bioactive scaffolds
           with scaffold networks: Improved compound set enrichment from primary screening data.
           Journal of Chemical Information and Modeling, 51(7), 1528–1538.
    .. [2] Varin, T., Gubler, H., Parker, C., Zhang, J., Raman, P., Ertl, P. and Schuffenhauer, A. (2010)
           Compound Set Enrichment: A Novel Approach to Analysis of Primary HTS Data.
           Journal of Chemical Information and Modeling, 50(12), 2067-2078.

    """
    set_node_attributes(scaffoldgraph, calc_scaffold_enrichment(scaffoldgraph, activity, mode, alternative, p))
    bonferroni = bonferroni_correction(scaffoldgraph, crit)
    result = []
    for scaffold, data in scaffoldgraph.get_scaffold_nodes(True):
        if data['pval'] < bonferroni[data['hierarchy']]:
            result.append((scaffold, data))
    return tuple(sorted(result, key=lambda x: x[1]['pval']))
