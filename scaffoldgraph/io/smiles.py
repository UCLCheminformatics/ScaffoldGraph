"""
scaffoldgraph.io.smiles
"""

from rdkit.Chem import SmilesMolSupplier

from .supplier import EnumeratedMolSupplier, MolSupplier


def read_smiles_file(smiles_file, delimiter=' ', smiles_column=0,
                     name_column=1, header=False, requires_length=False):
    """Read a SMILES file.

    Parameters
    ----------
    smiles_file: path to a SMILES file
    requires_length: If True returns an enumerated Mol
        supplier, i.e. when monitoring progress

    Returns
    -------
    either a MolSupplier or an EnumeratedSupplier
    depending on whether a length is required
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

    supplier = SmilesMolSupplier(smiles_file, delimiter, smiles_column, name_column, header, True)

    return EnumeratedMolSupplier(supplier, count)


def smiles_count(smiles_file):
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
