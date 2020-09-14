"""
scaffoldgraph.prioritization.rule_io

Functions for reading prioritization rulesets from a file.
Can be used to specify custom rulesets when using the CLI.
"""

from scaffoldgraph.prioritization.original_rules import *
from scaffoldgraph.prioritization.generic_rules import *


rule_name_to_class = {
    'scpabsdelta': SCPAbsDelta,
    'scpdelta': SCPDelta,
    'scpnumlinkerbonds': SCPNumLinkerBonds,
    'scpnumaromaticrings': SCPNumAromaticRings,
    'scpnumhetatoms': SCPNumHetAtoms,
    'scpnumnatoms': SCPNumNAtoms,
    'scpnumoatoms': SCPNumOAtoms,
    'scpnumsatoms': SCPNumSAtoms,
    'scpnumxatoms': SCPNumXAtoms,
    'rrpringsize': RRPRingSize,
    'rrplinkerlength': RRPLinkerLength,
    'rrphetatomlinked': RRPHetAtomLinked,
    'rrplinkerlengthx': RRPLinkerLengthX,
    'rrpnumhetatoms': RRPNumHetAtoms,
    'rrpnumnatoms': RRPNumNAtoms,
    'rrpnumoatoms': RRPNumOAtoms,
    'rrpnumsatoms': RRPNumSAtoms,
    'rrpnumxatoms': RRPNumXAtoms,
    'rrpringsizex': RRPRingSizeX,
    'rspabsdelta': RSPAbsDelta,
    'rspdelta': RSPDelta,
    'rspnumaromaticrings': RSPNumAromaticRings,
    'rspnumhetatoms': RSPNumHetAtoms,
    'rspnumnatoms': RSPNumNAtoms,
    'rspnumoatoms': RSPNumOAtoms,
    'rspnumsatoms': RSPNumSAtoms,
    'rspnumxatoms': RSPNumXAtoms,
    'tiebreaker': Tiebreaker,
    'originalrule01': OriginalRule01,
    'originalrule02': OriginalRule02,
    'originalrule03': OriginalRule03,
    'originalrule04': OriginalRule04,
    'originalrule05': OriginalRule05,
    'originalrule06': OriginalRule06,
    'originalrule07': OriginalRule07,
    'originalrule08': OriginalRule08,
    'originalrule09a': OriginalRule09a,
    'originalrule09b': OriginalRule09b,
    'originalrule09c': OriginalRule09c,
    'originalrule10': OriginalRule10,
    'originalrule11': OriginalRule11,
    'originalrule12': OriginalRule12,
    'originalrule13': OriginalRule13,
}


def read_rule_file(filename):
    """Read rules from a file.

    Parameters
    ----------
    filename : str
        Name of rule file.

    Returns
    -------
    list
        list of rule objects.

    Raises
    ------
    ValueError
        Raised if any of the rules defined in the
        rule file are not implemented.

    """
    rules = []
    with open(filename, 'r') as f:
        for line in f.readlines():
            tokens = line.strip().split('_')
            if len(tokens) == 0:
                continue
            rule_name = tokens[0]
            rule_cls = rule_name_to_class.get(
                rule_name.lower(), None)
            if rule_cls is None:
                raise ValueError(f'Rule {rule_name} is not defined')
            if len(tokens) > 2:
                rule = rule_cls(tokens[1], *list(map(int, tokens[2:])))
            else:
                rule = rule_cls(*tokens[1:])
            rules.append(rule)
    return rules
