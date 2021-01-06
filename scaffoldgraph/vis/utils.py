"""
scaffoldgraph.vis.utils
"""

import matplotlib.pyplot as plt
import matplotlib as mpl

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


def embed_node_mol_images(graph, size=(350, 300), draw_options=None, skip_existing=True):
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
    skip_existing : bool
        Skip node if it contains an 'img' attribute.
        The default is True.

    """
    for node, data in graph.nodes(data=True):
        if skip_existing and data.get('img', None):
            continue
        elif data.get('type', None) == 'scaffold':
            data['img'] = smiles_to_image(node, size, draw_options)
        elif data.get('type', None) == 'molecule':
            data['img'] = smiles_to_image(data['smiles'], size, draw_options)
        else:
            data['img'] = ''


def remove_node_mol_images(graph):
    """Remove embeded images from a graph.

    Parameters
    ----------
    graph : ScaffoldGraph
        Input ScaffoldGraph

    """
    for node, data in graph.nodes(data=True):
        _ = data.pop('img', None)


def rgba_to_hex(scalar_mappable, value):
    """str: rgba to hex."""
    rgba = scalar_mappable.to_rgba(value)
    c_hex = mpl.colors.to_hex(rgba, keep_alpha=False)
    return c_hex


def cmap_to_scalar_mappable(cmap, vmin, vmax):
    """Convert matplotlib Colormap to a ScalarMappable.

    Parameters
    ----------
    cmap : matplotlib.colors.Colormap
    vmin : float
        Minimum value for normalization.
    vmax : float
        Maximum value for normalization.

    Returns
    -------
    matplolib.cm.ScalarMappable

    """
    cnorm = mpl.colors.Normalize(vmin, vmax)
    scalar = mpl.cm.ScalarMappable(norm=cnorm, cmap=cmap)
    return scalar


def color_nodes_by_attribute(graph, attribute, cmap, node_type, label='color'):
    """
    Add an attribute to nodes in a ScaffoldGraph containing a color hex code,
    calculated from a paticular node attribute and a matplotlib cmap. The
    operation is perfomred in-place.

    Can be used for adding colors to ScaffoldGraph visualizations.

    Parameters
    ----------
    graph : ScaffoldGraph
        Input ScaffoldGraph
    attribute : str
        Key for the attibute from which to calculate a color.
    cmap : str or matplotlib.colors.Colormap
        A matplotlib cmap or name of a cmap e.g. 'BuPu' for
        calculating a nodes colour.
    node_type : str
        The type of node to process e.g. 'scaffold' / 'molecule'
    label : str, optional
        The attribute label to use for storing the color.
        The default is 'color'.

    """
    # Cmap may be a string or a Colormap
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    else:
        if not issubclass(type(cmap), mpl.colors.Colormap):
            raise ValueError('cmap must be a string or a matplotlib Colormap')

    # Get attribute range.
    _, attr = zip(*graph._get_nodes_with_type(node_type, attribute, None))
    attr = list(filter(lambda x: x is not None, attr))
    attr = list(map(float, attr))
    vmin, vmax = min(attr), max(attr)

    # Assign colors to each node.
    scalar_mappable = cmap_to_scalar_mappable(cmap, vmin, vmax)
    for node, data in graph._get_nodes_with_type(node_type, True, None):
        attr_val = data.get(attribute, None)
        if not attr_val:
            color = '#EEEEEE'  # Set a neutral default.
        else:
            color = rgba_to_hex(scalar_mappable, attr_val)
        data[label] = color


def color_scaffold_nodes_by_attribute(graph, attribute, cmap, label='color'):
    """
    Add an attribute to scaffold nodes in a ScaffoldGraph containing a color hex code,
    calculated from a paticular scaffold node attribute and a matplotlib cmap. The
    operation is perfomred in-place.

    Can be used for adding colors to ScaffoldGraph visualizations.

    Parameters
    ----------
    graph : ScaffoldGraph
        Input ScaffoldGraph
    attribute : str
        Key for the attibute from which to calculate a color.
    cmap : str or matplotlib.colors.Colormap
        A matplotlib cmap or name of a cmap e.g. 'BuPu' for
        calculating a nodes colour.
    label : str, optional
        The attribute label to use for storing the color.
        The default is 'color'.

    See Also
    --------
    color_molecule_nodes_by_attribute

    """
    color_nodes_by_attribute(graph, attribute, cmap, 'scaffold', label)


def color_molecule_nodes_by_attribute(graph, attribute, cmap, label='color'):
    """
    Add an attribute to molecule nodes in a ScaffoldGraph containing a color hex code,
    calculated from a paticular molecule node attribute and a matplotlib cmap. The
    operation is perfomred in-place.

    Can be used for adding colors to ScaffoldGraph visualizations.

    Parameters
    ----------
    graph : ScaffoldGraph
        Input ScaffoldGraph
    attribute : str
        Key for the attibute from which to calculate a color.
    cmap : str or matplotlib.colors.Colormap
        A matplotlib cmap or name of a cmap e.g. 'BuPu' for
        calculating a nodes colour.
    label : str, optional
        The attribute label to use for storing the color.
        The default is 'color'.

    See Also
    --------
    color_scaffold_nodes_by_attribute

    """
    color_nodes_by_attribute(graph, attribute, cmap, 'molecule', label)


def add_root_node(graph):
    """Add a root node to a scaffoldgraph.

    Parameters
    ----------
    graph : ScaffoldGraph
        Graph to add root node.

    """
    graph.add_node('root', type='root', hierarchy=0)
    edges = [('root', s) for s, d in graph.in_degree if d == 0 and s != 'root']
    graph.add_edges_from(edges, type=2)


def remove_root_node(graph):
    """Remove a root node from a scaffoldgraph.

    Parameters
    ----------
    graph : Scaffoldgraph
        Graph from which to remove root node.

    """
    if 'root' in graph:
        graph.remove_node('root')
