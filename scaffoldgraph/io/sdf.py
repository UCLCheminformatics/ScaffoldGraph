"""
scaffoldgraph.io.sdf

Contains functions for reading and writing from/to SDF.
"""

from rdkit.Chem import ForwardSDMolSupplier, SDWriter, MolFromSmiles

from .supplier import MolSupplier, EnumeratedMolSupplier


def read_sdf(sdf_file, requires_length=False):
    """Read molecules from an SDF.

    Parameters
    ----------
    sdf_file : file-like object
        An open SDF.
    requires_length : bool, optional
        If True returns an enumerated MolSupplier,
        i.e. when monitoring progress. The default
        is False.

    Returns
    -------
    MolSupplier or EnumeratedSupplier

    """
    supplier = ForwardSDMolSupplier(sdf_file)
    if not requires_length:
        return MolSupplier(supplier)
    count = sdf_count(sdf_file)
    sdf_file.seek(0)
    return EnumeratedMolSupplier(supplier, count)


def write_sdf_file(scaffold_graph, output_file):
    """Write an SDF file from a ScaffoldGraph.

    All scaffolds in the scaffoldgraph are written to the
    SDF, while molecules are ignored. Scaffolds are sorted
    in ascending order according to their hierarchy level.

    The output follows the standard SDF specification with
    the added property fields:

        TITLE field: scaffold ID
        SUBSCAFFOLDS field: list of sub-scaffold IDs
        HIERARCHY field: hierarchy level of scaffold
        SMILES field: scaffold canonical SMILES

    Parameters
    ----------
    scaffold_graph : scaffoldgraph.core.ScaffoldGraph
        ScaffoldGraph to be written to an SDF.
    output_file : str
        Filepath to an output file.

    """
    N = scaffold_graph.num_scaffold_nodes
    sorted_scaffolds = sorted(scaffold_graph.get_scaffold_nodes(data=True), key=lambda x: x[1]['hierarchy'])
    mapping = dict(zip([s[0] for s in sorted_scaffolds], range(0, N)))
    writer = SDWriter(output_file)
    for scaffold, data in sorted_scaffolds:
        molecule = MolFromSmiles(scaffold)
        if molecule is not None:
            subscaffolds = list(scaffold_graph.predecessors(scaffold))
            molecule.SetProp('_Name', mapping[scaffold])
            molecule.SetIntProp('HIERARCHY', scaffold_graph.nodes[scaffold]['HIERARCHY'])
            molecule.SetProp('SMILES', scaffold)
            molecule.SetProp('SUBSCAFFOLDS', ', '.join([str(mapping[s]) for s in subscaffolds]))
            writer.write(molecule)
    writer.close()


def sdf_count(file_obj):
    """Count the number of molecules in an SDF file.

    Counts the number of times '$$$$' occurs at the start of lines
    in the file.

    Parameters
    ----------
    file_obj : file-like object

    Returns
    -------
    int
        The number of molecules in the file.

    """
    return sum(1 for line in file_obj if line[:4] == b'$$$$')
