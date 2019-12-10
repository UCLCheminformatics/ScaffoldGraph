"""
scaffoldgraph.io.supplier
"""

from loguru import logger


class MolSupplier(object):
    """A wrapper for rdkit Mol suppliers"""

    def __init__(self, supplier):
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
    """A wrapper for rdkit Mol suppliers, providing the number of mols in the supplier,
    for use with progress monitoring"""

    def __init__(self, supplier, length):
        """Initialize an enumerated supplier

        Parameters
        ----------
        supplier: an rdkit Mol Supplier
        length: number of Mols in the supplier
        """

        super(EnumeratedMolSupplier, self).__init__(supplier)
        self.n = length

    def __len__(self):
        return self.n
