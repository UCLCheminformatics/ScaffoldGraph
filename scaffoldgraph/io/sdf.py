"""
scaffoldgraph.io.sdf
"""

from rdkit.Chem import ForwardSDMolSupplier, SDWriter, MolFromSmiles

from .supplier import MolSupplier, EnumeratedMolSupplier


def read_sdf(sdf_file, requires_length=False):
    """Read an sdf file.

    Parameters
    ----------
    sdf_file: A file-like object
    requires_length: If True returns an enumerated Mol
        supplier, i.e. when monitoring progress

    Returns
    -------
    either a MolSupplier or an EnumeratedSupplier
    depending on whether a length is required
    """

    supplier = ForwardSDMolSupplier(sdf_file)
    if not requires_length:
        return MolSupplier(supplier)
    count = sdf_count(sdf_file)
    sdf_file.seek(0)
    return EnumeratedMolSupplier(supplier, count)


def write_sdf_file(scaffold_graph, output_file):
    """Write an SDF file from a scaffoldgraph

    Parameters
    ----------
    scaffold_graph (sg.ScaffoldGraph): graph to be converted
    output_file (str): path to output file
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
    file_obj: A file-like object

    Returns
    -------
    count: The number of molecules in the file (int)
    """
    return sum(1 for line in file_obj if line[:4] == b'$$$$')
