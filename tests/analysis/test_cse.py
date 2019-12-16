"""
scaffoldgraph tests.analysis.test_cse
"""

import random
import scaffoldgraph as sg

from scaffoldgraph.analysis import CSE
from .. import mock_sdf

def test_cse(sdf_file):
    network = sg.ScaffoldNetwork.from_sdf(sdf_file)
    mol_nodes = network.get_molecule_nodes()

    mapping = {n: random.random() for n in mol_nodes}
    cse = CSE(hypothesis_test='ks')
    cse.fit(network, mapping, progress=False)
    assert len(cse) == 8
    assert 'PVAL' in cse['c1ccccc1'].keys()
    assert 'KS' in cse['c1ccccc1'].keys()

    df = cse.to_dataframe()
    assert 'SCAFFOLD' in df.columns
    assert 'PVAL' in df.columns
    assert 'CRIT' in df.columns

    mapping = {n: random.getrandbits(1) for n in mol_nodes}
    cse = CSE(hypothesis_test='binom')
    cse.fit(network, mapping, progress=False)
    assert len(cse) == 8
    assert 'PVAL' in cse['c1ccccc1'].keys()
    assert 'CRIT' in df.columns

    assert repr(cse) == '<CSE at {}>'.format(hex(id(cse)))
