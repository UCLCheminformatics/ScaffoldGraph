"""
scaffoldgraph.io.smiles

Contains functions for reading molecules from SMILES files.
"""

from rdkit.Chem import SmilesMolSupplier

from .supplier import EnumeratedMolSupplier, MolSupplier


def read_smiles_file(smiles_file, delimiter=' ', smiles_column=0,
                     name_column=1, header=False, requires_length=False):

    """Read molecules from a SMILES file.

    Parameters
    ----------
    smiles_file : str
        File path to a SMILES file.
    delimiter : str, optional
        Delimiter used in SMILES file. The default is ' '.
    smiles_column : int, optional
        SMILES column index. The default is 0.
    name_column : int, optional
        Molecule name/ID column index. The default is 1.
    header : bool, optional
        Whether the SMILES file contains a header.
        The default is False.
    requires_length : bool, optional
        If True returns an enumerated Mol supplier, i.e. when
        monitoring progress. The default is False.

    Returns
    -------
    MolSupplier or EnumeratedSupplier

    """
    if requires_length is False:
        return MolSupplier(
            SmilesMolSupplier(
                smiles_file,
                delimiter,
                smiles_column,
                name_column,
                header,
                True))

    count = smiles_count(smiles_file)
    if header is True:
        count -= 1

    supplier = SmilesMolSupplier(
        smiles_file, delimiter, smiles_column, name_column, header, True
    )

    return EnumeratedMolSupplier(supplier, count)


def smiles_count(smiles_file):
    """int : Return the number of lines in a SMILES file."""
    f = open(smiles_file, 'rb')
    lines = 0
    buf_size = 1024 * 1024
    read_f = f.read
    buf = read_f(buf_size)
    while buf:
        lines += buf.count(b'\n')
        buf = read_f(buf_size)
    f.close()
    return lines
