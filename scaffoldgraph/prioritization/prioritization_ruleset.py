"""
scaffoldgraph.prioritization.prioritization_ruleset

Implements a ruleset for scaffold prioritization when constructing scaffold trees.
"""

from .prioritization_rules import BaseScaffoldFilterRule


class ScaffoldRuleSet(object):
    """
    Class defining a set of rules used for scaffold prioritization.

    Rules added to the rule set must subclass the BaseScaffoldFilterRule.

    """
    def __init__(self, rules=None, name=None):
        """
        Initialize a rule set with an iterable of rules and an
        optional name.

        Parameters
        ----------
        rules : iterable, optional
            An iterable of rules. The default is None.
        name : str, optional
            Name of rule set. The default is None.

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
        """list : Return rules as a list."""
        return self._rules

    def filter_scaffolds(self, child, parents):
        """Filter a set of parent scaffolds using the defined rules.

        Method is called internally by scaffold graph constructors.
        __call__ is an alias for this function.

        Parameters
        ----------
        child : scaffoldgraph.core.Scaffold
            Child scaffold.
        parents : list
            Parent scaffolds.

        Returns
        -------
        parent : scaffoldgraph.core.Scaffold
            The scaffold retained after filtering.

        Raises
        ------
        ValueError
            Raised if the ruleset contains no rules.
        ValueError
            Raised if the iterable of parent scaffolds
            is empty.
        ValueError
            Raised if more than one scaffold is left after
            all of the filter rules are evaluated. The RuleSet
            may require a tie-breaker rule.

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
        raise ValueError('Filter error, more than one remaining scaffold '
                         'after filter rules applied. Rule set may require '
                         'a tie-breaker rule')

    def add_rule(self, rule):
        """Appends a rule to the ruleset.

        Parameters
        ----------
        rule : BaseScaffoldFilterRule
            Scaffold filter rule with base class ``BaseScaffoldFilterRule``.

        """
        if self.check_valid_rule(rule):
            self._rules.append(rule)
        else:
            raise TypeError('rule must be a subclass of BaseScaffoldRule')

    def insert_rule(self, rule, index):
        """Inserts a rule into the ruleset at supplied index.

        Parameters
        ----------
        rule : BaseScaffoldFilterRule
            Scaffold filter rule with base class ``BaseScaffoldFilterRule``.
        index : int
            Position in list to insert rule.

        """
        if self.check_valid_rule(rule):
            self._rules.insert(index, rule)
        else:
            raise TypeError('rule must be a subclass of BaseScaffoldRule')

    def delete_rule(self, index):
        """Deletes a rule from the ruleset at supplied index.

        Parameters
        ----------
        index : int
            Position in list to delete rule.

        """
        self._rules.__delitem__(index)

    @classmethod
    def from_rule_file(cls, filename, name=None):
        """Create a scaffold rule set from a rule set file.

        A rule set file is a text file specifying the names of
        rules to include in the ruleset seperated by new lines.
        The rule names must belong to either the original set
        or the generic set. The name of the rule corresponds to
        the class name of the desired rule. i.e. for OriginalRule01
        the file should contain the string OriginalRule01 followed
        by a new line. When including generic rules, min or max
        can be specified by including min or max after the name
        seperated by an underscore. i.e. RRPNumHetAtoms_min.
        For Rules which contain further arguments, these can be
        appended to the name with underscores. i.e.
        RRPRingSizeX_max_6. In this case the rule will prioritize
        scaffolds where the removed rings size is equal to 6.

        Parameters
        ----------
        filename : str
            File name of the rule set file.
        name : str, optional
            Name to assign rule set.

        See Also
        --------
        scaffoldgraph.prioritization.original_rules
        scaffoldgraph.prioritization.generic_rules

        """
        from .rule_io import read_rule_file
        rules = read_rule_file(filename)
        return cls(rules, name)

    @staticmethod
    def check_valid_rule(rule):
        """bool : Returns True if rule is a valid scaffold filter rule."""
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
