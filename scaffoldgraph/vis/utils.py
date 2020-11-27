"""
scaffoldgraph.vis.utils
"""

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem

from loguru import logger
from urllib import parse


def _maybe_kekulize(mol):
    """Private: attempt to kekulize a molecule."""
    try:
        Chem.Kekulize(mol)
    except Chem.KekulizeException:
        smi = Chem.MolToSmiles(mol)
        logger.warning(f'Failed to kekulize mol: {smi}')
    return mol


def smiles_to_svg(smiles, size=(350, 300), draw_options=None):
    """Create an SVG string from a SMILES string.

    Parameters
    ----------
    smiles : str
        SMILES to create SVG image.
    size : tuple, optional
        Size of image, the default is (350, 300).
    draw_options : rdMolDraw2D.MolDrawOptions
        Options to pass to the drawer.

    Returns
    -------
    svg : str
        SVG text for molecule.

    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return ''
    mol = _maybe_kekulize(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(*size)
    if draw_options:
        drawer.SetDrawOptions(draw_options)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def smiles_to_image(smiles, size=(350, 300), draw_options=None):
    """Create an SVG image from a SMILES string (ready for HTML).

    Parameters
    ----------
    smiles : str
        SMILES to create SVG image.
    size : tuple, optional
        Size of image, the default is (350, 300).
    draw_options : rdMolDraw2D.MolDrawOptions
        Options to pass to the drawer.

    Returns
    -------
    svg : str
        SVG image path.

    """
    svg = smiles_to_svg(smiles, size, draw_options)
    img_path = 'data:image/svg+xml;charset=utf-8,'
    img_path += parse.quote(svg, safe='')
    return img_path


def embed_node_mol_images(graph, size=(350, 300), draw_options=None):
    """Embed molecule images into a graph.

    Images are added as an attribute 'img' to each node with an
    available SMILES string ('molecule', 'scaffold'). The graph
    is modified in-place.

    Parameters
    ----------
    graph : ScaffoldGraph
        Input ScaffoldGraph.
    size : tuple, optional
        Size of image, the default is (350, 300).
    draw_options : rdMolDraw2D.MolDrawOptions
        Options to pass to the drawer.

    """
    for node, data in graph.nodes(data=True):
        if data.get('type', None) == 'scaffold':
            data['img'] = smiles_to_image(node, size, draw_options)
        elif data.get('type', None) == 'molecule':
            data['img'] = smiles_to_image(data['smiles'], size, draw_options)
        else:
            data['img'] = ''
