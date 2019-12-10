"""
scaffoldgraph.prioritization.original_rules

Implements rules from the paper:
'The Scaffold Tree − Visualization of the Scaffold Universe by Hierarchical Scaffold Classification'
"""

from itertools import chain

from rdkit.Chem import MolFromSmarts

from .prioritization_rules import *
from .prioritization_ruleset import ScaffoldRuleSet


class OriginalRule01(ScaffoldFilterRule):
    """Remove heterocycles of size 3 first"""

    def condition(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        ring_atomic_nums = [a.GetAtomicNum() for a in removed_ring.atoms]
        ring_num_het = len([a for a in ring_atomic_nums if a != 1 and a != 6])
        return removed_ring.size == 3 and ring_num_het == 1

    @property
    def name(self):
        return 'original rule 01'


class OriginalRule02(ScaffoldFilterRule):
    """Do not remove rings with >= 12 atoms if there are smaller rings to remove"""

    def condition(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        return removed_ring.size < 12

    @property
    def name(self):
        return 'original rule 02'


class OriginalRule03(ScaffoldMinFilterRule):
    """Choose the parent scaffold with the smallest number of acyclic linker bonds"""

    acyc_linker_smarts = MolFromSmarts('*!@!=!#*')

    def get_property(self, child, parent):
        matches = parent.mol.GetSubstructMatches(self.acyc_linker_smarts)
        return len(matches)

    @property
    def name(self):
        return 'original rule 03'


class OriginalRule04(ScaffoldMaxFilterRule):
    """Retain bridged rings, spiro rings and nonlinear fusion patterns with preference"""

    def get_property(self, child, parent):
        nr = parent.rings.count
        rb = list(chain(*parent.rings.bond_rings))
        nrrb = len(rb) - len(set(rb))
        return abs(nrrb - (nr - 1))

    @property
    def name(self):
        return 'original rule 04'


class OriginalRule05(ScaffoldFilterRule):
    """Bridged ring systems retained with preference over spiro rings,
    Rings with a positive signed delta are retained"""

    def condition(self, child, parent):
        nr = parent.rings.count
        rb = list(chain(*parent.rings.bond_rings))
        nrrb = len(rb) - len(set(rb))
        delta = nrrb - (nr - 1)
        return delta >= 1

    @property
    def name(self):
        return 'original rule 05'


class OriginalRule06(ScaffoldFilterRule):
    """Remove rings of size 3, 5 and 6 first"""

    def condition(self, child, parent):
        rr_size = child.rings[parent.removed_ring_idx].size
        return rr_size == 3 or rr_size == 5 or rr_size == 6

    @property
    def name(self):
        return 'original rule 06'


class OriginalRule07(BaseScaffoldFilterRule):
    """A fully aromatic ring system must not be dissected in a way that the resulting system
     is not aromatic anymore.

    UNIMPLEMENTED
    This is tricky to implement and should probably be done during the fragmentation process,
    although for efficiency it might be better to ignore this rule, rdkit seems to catch many
    of these cases in the partial sanitization as we do not attempt to change atom types when
    this event occurs. (SNG also skips this step)
    """

    def filter(self, child, parents):
        return parents

    @property
    def name(self):
        return 'original rule 07'


class OriginalRule08(ScaffoldMinFilterRule):
    """Remove rings with the least hetero atoms first"""

    def get_property(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        ring_atomic_nums = [a.GetAtomicNum() for a in removed_ring.atoms]
        return len([a for a in ring_atomic_nums if a != 1 and a != 6])

    @property
    def name(self):
        return 'original rule 08'


class OriginalRule09a(ScaffoldMinFilterRule):
    """Remove scaffolds with least nitrogen atoms in deleted ring"""

    def get_property(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        ring_atomic_nums = [a.GetAtomicNum() for a in removed_ring.atoms]
        return ring_atomic_nums.count(7)

    @property
    def name(self):
        return 'original rule 09a'


class OriginalRule09b(ScaffoldMinFilterRule):
    """Remove scaffolds with least oxygen atoms in deleted ring"""

    def get_property(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        ring_atomic_nums = [a.GetAtomicNum() for a in removed_ring.atoms]
        return ring_atomic_nums.count(8)

    @property
    def name(self):
        return 'original rule 09b'


class OriginalRule09c(ScaffoldMinFilterRule):
    """Remove scaffolds with least sulphur atoms in deleted ring"""

    def get_property(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        ring_atomic_nums = [a.GetAtomicNum() for a in removed_ring.atoms]
        return ring_atomic_nums.count(16)

    @property
    def name(self):
        return 'original rule 09c'


class OriginalRule10(ScaffoldMinFilterRule):
    """Smaller rings are removed first"""

    def get_property(self, child, parent):
        return child.rings[parent.removed_ring_idx].size

    @property
    def name(self):
        return 'original rule 10'


class OriginalRule11(ScaffoldFilterRule):
    """Retain non-aromatic rings with preference"""

    def condition(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        return all([bond.GetIsAromatic() for bond in removed_ring.bonds])

    @property
    def name(self):
        return 'original rule 11'


class OriginalRule12(ScaffoldFilterRule):
    """Remove rings first where the linker is attached to a ring hetero atom at either end of the linker"""

    connection_patt = MolFromSmarts('[R]!@!=*')

    def condition(self, child, parent):
        removed_ring = child.rings[parent.removed_ring_idx]
        connections = {x[0] for x in child.mol.GetSubstructMatches(self.connection_patt)}
        ring_connections = connections.intersection(removed_ring.aix)
        connection_atomic_nums = [child.atoms[x].GetAtomicNum() for x in ring_connections]
        return len([a for a in connection_atomic_nums if a != 1 and a != 6]) > 0

    @property
    def name(self):
        return 'original rule 12'


class OriginalRule13(BaseScaffoldFilterRule):
    """Tie-breaker rule (alphabetical)"""

    def filter(self, child, parents):
        return [sorted(parents, key=lambda p: p.smiles)[0]]

    @property
    def name(self):
        return 'original rule 13'


all_rules = [
    OriginalRule01(),
    OriginalRule02(),
    OriginalRule03(),
    OriginalRule04(),
    OriginalRule05(),
    OriginalRule06(),
    OriginalRule07(),
    OriginalRule08(),
    OriginalRule09a(),
    OriginalRule09b(),
    OriginalRule09c(),
    OriginalRule10(),
    OriginalRule11(),
    OriginalRule12(),
    OriginalRule13(),
]

original_ruleset = ScaffoldRuleSet(all_rules, name='Original Rules')
