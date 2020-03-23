"""
scaffoldgraph.utils.misc

Defines miscellaneous functions used within scaffoldgraph
"""

from rdkit import Chem


def canonize_smiles(smiles, failsafe=True):
    """Canonize a SMILES string (with failsafe)"""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None and failsafe:
        return smiles
    return Chem.MolToSmiles(mol)
