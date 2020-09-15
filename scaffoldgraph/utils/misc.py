"""
scaffoldgraph.utils.misc

Defines miscellaneous functions used within scaffoldgraph.
"""

from rdkit import Chem


def canonize_smiles(smiles, failsafe=True):
    """Canonize a SMILES string (with failsafe).

    Parameters
    ----------
    smiles : str
        SMILES string to canonize.
    failsafe : bool
        If True, if the SMILES fails to parse
        the input SMILES is returned instead
        of raising an error.

    Returns
    -------
    str
        The canonical SMILES representation.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None and failsafe:
        return smiles
    return Chem.MolToSmiles(mol)
