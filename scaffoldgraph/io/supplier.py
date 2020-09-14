"""
scaffoldgraph.io.supplier

Contains utilities for io within scaffoldgraph.
"""

from loguru import logger


class MolSupplier(object):
    """A wrapper for rdkit Mol suppliers.

    Provides logging for molecule parsing errors in
    a way that is compatible with scaffoldgraphs
    logging system.

    Notes
    -----
    Technically the supplier can be used with any iterable python object
    containing or supplying rdkit Mol objects.

    See Also
    --------
    EnumeratedMolSupplier

    """
    def __init__(self, supplier):
        """Initialize an EnumeratedMolSupplier.

        Parameters
        ----------
        supplier : iterable
            An rdkit Mol Supplier.

        """
        self.supplier = supplier
        self.cursor = 1

    def __iter__(self):
        return self

    def __next__(self):
        mol = next(self.supplier)
        if mol is None:
            logger.warning('Molecule {} could not be parsed'.format(
                self.cursor
            ))
        self.cursor += 1
        return mol


class EnumeratedMolSupplier(MolSupplier):
    """
    A wrapper for rdkit Mol suppliers, providing the number of mols in the supplier,
    for use with progress monitoring.

    Attributes
    ----------
    n : int
        The length of the supplier

    Notes
    -----
    Technically the supplier can be used with any iterable python object
    containing or supplying rdkit Mol objects.

    See Also
    --------
    MolSupplier

    """
    def __init__(self, supplier, length):
        """Initialize an EnumeratedMolSupplier.

        Parameters
        ----------
        supplier : iterable
            An rdkit Mol Supplier.
        length : int
            Number of Mols in the supplier.

        """
        super(EnumeratedMolSupplier, self).__init__(supplier)
        self.n = length

    def __len__(self):
        return self.n
