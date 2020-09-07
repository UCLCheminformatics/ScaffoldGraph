"""
scaffoldgraph.prioritization.prioritization_ruleset

Implements a ruleset for scaffold prioritization when constructing scaffold trees
"""

from .prioritization_rules import BaseScaffoldFilterRule


class ScaffoldRuleSet(object):
    """
    Class defining a set of rules used for scaffold prioritization
    Rules added to the rule set must subclass the BaseScaffoldFilterRule

    Parameters
    ----------
    rules: an iterable of rules (optional, default: None)
    name: name of rule set (optional, default: 'ScaffoldRuleSet')
    """

    def __init__(self, rules=None, name=None):

        """Initialize a rule set with an iterable of rules and a name

        Parameters
        ----------
        rules: an iterable of rules (optional, default: None)
        name: name of rule set (optional, default: 'ScaffoldRuleSet')
        """

        self._rules = []
        if rules is not None:
            for rule in rules:
                self.add_rule(rule)
        self.name = name if name else 'ScaffoldRuleSet'

    def __call__(self, child, parents):
        return self.filter_scaffolds(child, parents)

    @property
    def rules(self):
        """Return rules as a list"""
        return self._rules

    def filter_scaffolds(self, child, parents):
        """Filter a set of parent scaffolds using the defined rules

        Method is called internally by scaffold graph constructors

        Parameters
        ----------
        child (Scaffold) : child scaffold
        parents (list) : list of parent scaffolds
        """

        if len(self) == 0:
            raise ValueError('No rules defined in rule set')
        if len(parents) == 0:
            raise ValueError('No parent scaffolds supplied to filter')
        elif len(parents) == 1:
            parent = parents.pop()
            parent.prioritization_rule = 'last remaining'
            return parent
        remaining = list(parents)
        for rule in self:
            filtered = rule.filter(child, remaining)
            if filtered:
                remaining = filtered
            if len(remaining) == 1:
                parent = remaining.pop()
                parent.prioritization_rule = rule.name
                return parent
        raise RuntimeError('Filter error, more than one remaining scaffold '
                           'after filter rules applied. Rule set may require '
                           'a tie-breaker rule')

    def add_rule(self, rule):
        """Appends a rule to the ruleset

        Parameters
        ----------
        rule : scaffold filter rule with base class BaseScaffoldFilterRule
        """

        if self.check_valid_rule(rule):
            self._rules.append(rule)
        else:
            raise TypeError('rule must be a subclass of BaseScaffoldRule')

    def insert_rule(self, rule, index):
        """Inserts a rule into the ruleset at supplied index

        Parameters
        ----------
        rule : scaffold filter rule with base class BaseScaffoldFilterRule
        index (int) : position in list to insert rule
        """
        if self.check_valid_rule(rule):
            self._rules.insert(index, rule)
        else:
            raise TypeError('rule must be a subclass of BaseScaffoldRule')

    def delete_rule(self, index):
        """Deletes a rule from the ruleset at supplied index

        Parameters
        ----------
        index (int) : position in list to delete rule
        """
        self._rules.__delitem__(index)

    @classmethod
    def from_rule_file(cls, filename, name=None):
        """Create a scaffold rule set from a rule set file

        Parameters
        ----------
        filename (str) : file name of the rule set file
        name (str, optional) : name to assign rule set
        """
        from .rule_io import read_rule_file
        rules = read_rule_file(filename)
        return cls(rules, name)

    @staticmethod
    def check_valid_rule(rule):
        """Returns True if rule is a valid scaffold filter rule"""
        return BaseScaffoldFilterRule in rule.__class__.__mro__

    def __getitem__(self, index):
        return self._rules[index]

    def __setitem__(self, index, rule):
        if self.check_valid_rule(rule):
            self._rules.__setitem__(index, rule)
        raise TypeError('rule must be a subclass of BaseScaffoldRule')

    def __len__(self):
        return len(self._rules)

    def __repr__(self):
        return '<{_cls} at {address}>'.format(
            _cls=self.__class__.__name__,
            address=hex(id(self))
        )
