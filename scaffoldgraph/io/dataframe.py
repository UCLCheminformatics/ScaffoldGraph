"""
scaffoldgraph.io.dataframe

Contains functions for reading molecules from pandas dataframes.
"""

from rdkit.Chem import MolFromSmiles
from loguru import logger


class DataFrameMolSupplier(object):
    """Class supplying rdkit Mols from a pandas DataFrame."""

    def __init__(self, df, smiles_column, name_column, data_cols=None):
        """Initialize DataFrameMolSupplier.

        Parameters
        ----------
        df : pandas.DataFrame
            Dataframe to read molecules from.
        smiles_column : str
            Key of column containing SMILES strings.
        name_column : str
            Key of column containing molecule name strings.
        data_cols : list, optional
            A list of column keys containg data to retain
            in molecule graph nodes. The default is None.

        """
        self.data_cols = data_cols
        if data_cols is None:
            self.supplier = zip(
                df[smiles_column].values,
                df[name_column].values
            )
        else:
            self.supplier = zip(
                df[smiles_column].values,
                df[name_column].values,
                df[data_cols].values
            )
        self.n = len(df[smiles_column])
        self.cursor = 1

    def __iter__(self):
        return self

    def __next__(self):
        values = next(self.supplier)
        try:
            mol = MolFromSmiles(values[0])
            mol.SetProp('_Name', str(values[1]))
            if self.data_cols is not None:
                for key, value in zip(self.data_cols, values[2]):
                    mol.SetProp(str(key), str(value))
        except AttributeError:
            logger.warning('Molecule {} : {} could not be parsed'.format(
                self.cursor, values[0]
            ))
            self.cursor += 1
            return None

        self.cursor += 1
        return mol

    def __len__(self):
        return self.n


def read_dataframe(df, smiles_column, name_column, data_columns=None):
    """Read molecules from a dataframe.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe to read molecules from.
    smiles_column : str
        Key of column containing SMILES strings.
    name_column : str
        Key of column containing molecule name strings.
    data_columns : list, optional
        A list of column keys containg data to retain
        in molecule graph nodes. The default is None.

    Returns
    -------
    DataFrameMolSupplier

    """
    return DataFrameMolSupplier(df, smiles_column, name_column, data_columns)
