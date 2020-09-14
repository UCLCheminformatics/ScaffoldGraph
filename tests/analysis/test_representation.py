"""
scaffoldgraph tests.analysis.test_representation
"""

from scaffoldgraph.analysis import calc_average_pairwise_similarity, get_over_represented_scaffold_classes
from ..test_network import long_test_network


def test_representation(network):
    aps = calc_average_pairwise_similarity(network)
    entry = list(aps.items())[0]
    assert entry[0] in network
    assert 'members' in entry[1]
    assert 'aps' in entry[1]
    assert type(entry[1]['aps']) == float
    assert type(entry[1]['members']) == int
    over = get_over_represented_scaffold_classes(network, 0.80)
    for entry in over:
        assert entry[1]['aps'] >= 0.80
    over = get_over_represented_scaffold_classes(network, 0.75, skip_aps=True)
    for entry in over:
        assert entry[1]['aps'] >= 0.75
