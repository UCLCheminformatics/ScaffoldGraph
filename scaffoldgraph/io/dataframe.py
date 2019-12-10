"""
scaffoldgraph.io.dataframe
"""

from rdkit.Chem import MolFromSmiles


class DataFrameMolSupplier(object):
    """"""

    def __init__(self, df, smiles_column, name_column):
        self.supplier = zip(df[smiles_column].values, df[name_column].values)
        self.n = len(df[smiles_column])
        self.cursor = 1

    def __iter__(self):
        return self

    def __next__(self):
        smiles, name = next(self.supplier)

        try:
            mol = MolFromSmiles(smiles)
            mol.SetProp('_Name', str(name))

        except AttributeError:
            logger.warning('Molecule {} could not be parsed'.format(
                self.cursor
            ))
            self.cursor += 1
            return None

        self.cursor += 1
        return mol

    def __len__(self):
        return self.n


def read_dataframe(df, smiles_column, name_column):
    """Read molecules from a dataframe.

    Parameters
    ----------

    Returns
    -------
    DataFrameMolSupplier
    """
    return DataFrameMolSupplier(df, smiles_column, name_column)
