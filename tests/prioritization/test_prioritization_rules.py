"""
scaffoldgraph tests.prioritization.test_prioritization_rules
"""

import pytest

from scaffoldgraph.prioritization.prioritization_rules import *


class MockScaffoldFilterRule(BaseScaffoldFilterRule):
    def filter(self, child, parents):
        return parents[1:]

    @property
    def name(self):
        return 'mock'


def test_prioritization_rules():
    """Test abstract ruletypes cannot be initialized."""

    with pytest.raises(TypeError):
        BaseScaffoldFilterRule()
    with pytest.raises(TypeError):
        ScaffoldFilterRule()
    with pytest.raises(TypeError):
        ScaffoldMinFilterRule()
    with pytest.raises(TypeError):
        ScaffoldMaxFilterRule()


def test_base_rule_subclass():
    """Test base class can be subclassed"""

    mock = MockScaffoldFilterRule()
    parents = [0, 1, 2, 3, 4]
    assert mock.name == 'mock'
    assert str(mock) == 'mock'
    assert mock.filter == [1, 2, 3, 4]
    assert mock(parents) == mock.filter(parents)
    assert repr(mock) == '<MockScaffoldFilterRule at {}>'.format(hex(id(mock)))
