"""
scaffoldgraph tests.prioritization.test_prioritization_ruleset
"""

import pytest

from scaffoldgraph.prioritization import *


@pytest.fixture(name='null_set')
def empty_ruleset():
    return ScaffoldRuleSet()


@pytest.fixture(name='original')
def original_ruleset():
    return original_ruleset


def test_empty_filter(null_set, scaffolds):
    """Test filtering with an empty ruleset raises an error"""
    with pytest.raises(ValueError):
        assert null_set.filter(null_set, scaffolds)


def test_name(null_set):
    assert null_set.name is None
    null_set.name = 'some_name'
    assert null_set.name == 'some_name'
