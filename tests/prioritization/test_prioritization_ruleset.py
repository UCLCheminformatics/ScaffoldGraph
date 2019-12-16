"""
scaffoldgraph tests.prioritization.test_prioritization_ruleset
"""

import pytest

from scaffoldgraph.prioritization.original_rules import *


@pytest.fixture(name='null_set')
def empty_ruleset():
    return ScaffoldRuleSet()


def test_empty_filter(null_set):
    """Test filtering with an empty ruleset raises an error"""
    with pytest.raises(ValueError):
        assert null_set.filter_scaffolds(null_set, [])
    with pytest.raises(ValueError):
        assert original_ruleset.filter_scaffolds('', [])


def test_name(null_set):
    assert null_set.name is 'ScaffoldRuleSet'
    null_set.name = 'some_name'
    assert null_set.name == 'some_name'


def test_rules():
    rules = original_ruleset.rules
    assert all([issubclass(x.__class__, BaseScaffoldFilterRule) for x in rules])


def test_builtins():
    assert len(original_ruleset) == 15
    assert isinstance(original_ruleset[0], BaseScaffoldFilterRule)
    assert original_ruleset.check_valid_rule(OriginalRule10())
    original_ruleset.add_rule(OriginalRule10())
    assert len(original_ruleset) == 16
    original_ruleset.insert_rule(OriginalRule10(), 1)
    assert original_ruleset[1].__class__ == OriginalRule10
    original_ruleset.delete_rule(16)
    original_ruleset.delete_rule(1)
    assert len(original_ruleset) == 15
    assert repr(original_ruleset) == '<ScaffoldRuleSet at {}>'.format(hex(id(original_ruleset)))


def test_errors(null_set):
    with pytest.raises(TypeError):
        null_set.add_rule('')
        null_set.insert_rule('', 0)
        null_set[0] = ''
