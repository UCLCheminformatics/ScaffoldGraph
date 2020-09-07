"""
scaffoldgraph tests.prioritization.test_generic_rules
"""

from rdkit import Chem

from scaffoldgraph.prioritization.generic_rules import *
from scaffoldgraph import get_next_murcko_fragments
from scaffoldgraph.utils import canonize_smiles
from scaffoldgraph.core import Scaffold


def fragment_and_filter(mol, rule):
    scf = Scaffold(mol)
    parents = get_next_murcko_fragments(mol)
    parents = list(map(lambda x: Scaffold(x), parents))
    output = [x.smiles for x in rule.filter(scf, parents)]
    return output


def _test_rule(mol, rule, expected):
    result = fragment_and_filter(mol, rule)
    assert len(result) == 1
    assert result[0] == expected


def _test_rule_min_max(mol, rule, expected_min, expected_max):
    _test_rule(mol, rule('min'), expected_min)
    _test_rule(mol, rule('max'), expected_max)


"""
SCP RULES scaffold property (parent)
------------------------------------
SCPNumLinkerBonds
SCPDelta
SCPAbsDelta
SCPNumAromaticRings
SCPNumHetAtoms
SCPNumNAtoms
SCPNumOAtoms
SCPNumSAtoms
"""


def test_scp_num_linker_bonds():
    test_smiles = 'O=C(NCCCCN1CCN(c2ccccc2)CC1)c1ccc2c(c1)Cc1ccccc1-2'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(NCCCCN1CCNCC1)c1ccc2c(c1)Cc1ccccc1-2')
    max_result = canonize_smiles('O=C(NCCCCN1CCN(c2ccccc2)CC1)c1ccc2c(c1)CC=C2')
    _test_rule_min_max(test_mol, SCPNumLinkerBonds, min_result, max_result)


def test_scp_delta():
    test_smiles = 'OC5CC31CC5CC1C4CCc2occc2C4CC3'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('OC1CCC2(CCC3c4ccoc4CCC3C2)C1')  # retain spiro
    max_result = canonize_smiles('OC1CC23CCC4C=CCCC4C2CC1C3')  # retain bridged
    _test_rule_min_max(test_mol, SCPDelta, min_result, max_result)


def test_scp_abs_delta():
    test_smiles = 'C1CC2CN3C(CC=CC3=O)C4C2N(C1)CCC4'
    test_mol = Chem.MolFromSmiles(test_smiles)
    max_result = canonize_smiles('C1CC2CNCC3CCCN(C1)C23')
    _test_rule(test_mol, SCPAbsDelta('max'), max_result)


def test_scp_num_het_atoms():
    test_smiles = 'C2Oc3ccccc3C(=O)C2=O'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('c1ccccc1')
    max_result = canonize_smiles('O=C1C=COCC1=O')
    _test_rule_min_max(test_mol, SCPNumHetAtoms, min_result, max_result)


def test_scp_num_aromatic_rings():
    test_smiles = 'N1CCCC(OC(=O)C(O)(c2ccccc2)c2ccccc2)C1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(OC1CCCNC1)C(O)c1ccccc1')
    max_result = canonize_smiles('OC(c1ccccc1)c1ccccc1')
    _test_rule_min_max(test_mol, SCPNumAromaticRings, min_result, max_result)


def test_scp_num_oxygen_atoms():
    test_smiles = 'C2Oc3ccccc3C(=O)C2=O'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('c1ccccc1')
    max_result = canonize_smiles('O=C1C=COCC1=O')
    _test_rule_min_max(test_mol, SCPNumOAtoms, min_result, max_result)


def test_scp_num_nitrogen_atoms():
    test_smiles = 'N1CCCC(OC(=O)C(O)(c2ccccc2)c2ccccc2)C1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('OC(c1ccccc1)c1ccccc1')
    max_result = canonize_smiles('O=C(OC1CCCNC1)C(O)c1ccccc1')
    _test_rule_min_max(test_mol, SCPNumNAtoms, min_result, max_result)


def test_scp_num_sulphur_atoms():
    test_smiles = 'c1csc2c(NCN3CCN(CCc4ccccc4)CC3)ncnc12'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('c1ccc(CCN2CCN(CNc3ccncn3)CC2)cc1')
    max_result = canonize_smiles('c1nc(NCN2CCNCC2)c2sccc2n1')
    _test_rule_min_max(test_mol, SCPNumSAtoms, min_result, max_result)


"""
RRP RULES removed ring property
-------------------------------
RRPRingSize
RRPLinkerLength
RRPHetAtomLinked
RRPNumHetAtoms
RRPNumNAtoms
RRPNumOAtoms
RRPNumSAtoms
"""


def test_rrp_ring_size():
    test_smiles = 'n1nc(-c2ccccc2)nc1=S'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('c1ccccc1')
    max_result = canonize_smiles('S=C1N=CN=N1')
    _test_rule_min_max(test_mol, RRPRingSize, min_result, max_result)


def test_rrp_linker_length():
    test_smiles = 'O=C1c2ccccc2-c2c(NCCc3ccccc3)c(=O)[nH]c3cccc1c23'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C1C=Cc2c(NCCc3ccccc3)c(=O)[nH]c3cccc1c23')
    max_result = canonize_smiles('O=C1c2ccccc2-c2cc(=O)[nH]c3cccc1c23')
    _test_rule_min_max(test_mol, RRPLinkerLength, min_result, max_result)


def test_rrp_het_atom_linked():
    test_smiles = 'O=C(NCc1ccccc1)N1CCN2C(=O)OC(c3ccccc3)(c3ccccc3)[C@@H]2C1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(NCc1ccccc1)N1CCN2C(=O)OC(c3ccccc3)[C@@H]2C1')
    max_result = canonize_smiles('O=C1OC(c2ccccc2)(c2ccccc2)[C@@H]2CNCCN12')
    _test_rule_min_max(test_mol, RRPHetAtomLinked, min_result, max_result)


def test_rrp_num_het_atoms():
    test_smiles = 'c1cccc2c(=O)[nH][nH]c(=O)c12'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=c1ccc(=O)[nH][nH]1')
    max_result = canonize_smiles('c1ccccc1')
    _test_rule_min_max(test_mol, RRPNumHetAtoms, min_result, max_result)


def test_rrp_num_nitrogen_atoms():
    test_smiles = 'c1cccc2c(=O)[nH][nH]c(=O)c12'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=c1ccc(=O)[nH][nH]1')
    max_result = canonize_smiles('c1ccccc1')
    _test_rule_min_max(test_mol, RRPNumNAtoms, min_result, max_result)


def test_rrp_num_oxygen_atoms():
    test_smiles = 'C1OC(=O)C2=C1CCC=C2'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C1C=CCO1')
    max_result = canonize_smiles('C1=CCCC=C1')
    _test_rule_min_max(test_mol, RRPNumOAtoms, min_result, max_result)


def test_rrp_num_sulphur_atoms():
    test_smiles = 'C1CSC(=NNC(=O)C(=O)CC2CCOCC2)N1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('N=C1NCCS1')
    max_result = canonize_smiles('C1CCOCC1')
    _test_rule_min_max(test_mol, RRPNumSAtoms, min_result, max_result)


"""
RSP Rules property of the ring system of a removed ring before removal
----------------------------------------------------------------------
RSPAbsDelta
RSPDelta
RSPNumAromaticRings
RSPNumHetAtoms
RSPNumNAtoms
RSPNumOAtoms
RSPNumRings
RSPNumSAtoms
"""


def test_rsp_delta():
    test_smiles = 'O=C1N(CCCC3CCNCC3)CCC12CCN1CCCC12'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C1N(CCCC2CCNCC2)CCC12CCNC2')
    max_result = canonize_smiles('O=C1NCCC12CCN1CCCC12')
    _test_rule_min_max(test_mol, RSPDelta, min_result, max_result)


def test_rsp_abs_delta():
    test_smiles = 'O=C1N(CCCC3CCNCC3)CCC12CCN1CCCC12'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C1NCCC12CCN1CCCC12')
    max_result = canonize_smiles('O=C1N(CCCC2CCNCC2)CCC12CCNC2')
    _test_rule_min_max(test_mol, RSPAbsDelta, min_result, max_result)


def test_rsp_num_aromatic_rings():
    test_smiles = 'O=C(c1c2ccccc2cc2ccccc12)N1CCC(N2CCC[C@@H](C(=O)N3CCOCC3)C2)CC1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(c1c2ccccc2cc2ccccc12)N1CCC(N2CCCCC2)CC1')
    max_result = canonize_smiles('O=C(c1cccc2ccccc12)N1CCC(N2CCC[C@@H](C(=O)N3CCOCC3)C2)CC1')
    _test_rule_min_max(test_mol, RSPNumAromaticRings, min_result, max_result)


def test_rsp_num_het_atoms():
    test_smiles = 'c1nc2ccc3nc(NC(=O)C(c4ccccc4)c4ccccc4)sc3c2s1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(Cc1ccccc1)Nc1nc2ccc3ncsc3c2s1')
    max_result = canonize_smiles('O=C(Nc1nc2ccccc2s1)C(c1ccccc1)c1ccccc1')
    _test_rule_min_max(test_mol, RSPNumHetAtoms, min_result, max_result)


def test_rsp_num_nitogen_atoms():
    test_smiles = 'c1nc2ccc3nc(NC(=O)C(c4ccccc4)c4ccccc4)sc3c2s1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(Cc1ccccc1)Nc1nc2ccc3ncsc3c2s1')
    max_result = canonize_smiles('O=C(Nc1nc2ccccc2s1)C(c1ccccc1)c1ccccc1')
    _test_rule_min_max(test_mol, RSPNumNAtoms, min_result, max_result)


def test_rsp_num_oxygen_atoms():
    test_smiles = 'c1c2c(c3occ(-c4ccccc4)c(=O)c3c1)C=CCO2'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=c1ccoc2c3c(ccc12)OCC=C3')
    max_result = canonize_smiles('O=c1c(-c2ccccc2)coc2ccccc12')
    _test_rule_min_max(test_mol, RSPNumOAtoms, min_result, max_result)


def test_rsp_num_sulphur_atoms():
    test_smiles = 'c1nc2ccc3nc(NC(=O)C(c4ccccc4)c4ccccc4)sc3c2s1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(Cc1ccccc1)Nc1nc2ccc3ncsc3c2s1')
    max_result = canonize_smiles('O=C(Nc1nc2ccccc2s1)C(c1ccccc1)c1ccccc1')
    _test_rule_min_max(test_mol, RSPNumSAtoms, min_result, max_result)


def test_rsp_num_rings():
    test_smiles = 'O=C(c1c2ccccc2cc2ccccc12)N1CCC(N2CCC[C@@H](C(=O)N3CCOCC3)C2)CC1'
    test_mol = Chem.MolFromSmiles(test_smiles)
    min_result = canonize_smiles('O=C(c1c2ccccc2cc2ccccc12)N1CCC(N2CCCCC2)CC1')
    max_result = canonize_smiles('O=C(c1cccc2ccccc12)N1CCC(N2CCC[C@@H](C(=O)N3CCOCC3)C2)CC1')
    _test_rule_min_max(test_mol, RSPNumRings, min_result, max_result)
