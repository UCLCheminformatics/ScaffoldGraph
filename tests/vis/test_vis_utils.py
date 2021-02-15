"""
scaffoldgraph tests.vis.test_vis_utils
"""

import scaffoldgraph.vis.utils as vis_utils
import matplotlib.pyplot as plt
import random
import pytest
import re

from rdkit.Chem.Draw import rdMolDraw2D
from rdkit import Chem

from scaffoldgraph.utils import suppress_rdlogger
from . import long_test_network


SVG_PATTERN = r'(?:<\?xml\b[^>]*>[^<]*)?(?:<!--.*?-->[^<]*)*(?:<svg|<!DOCTYPE svg)\b'
SVG_REGEX = re.compile(SVG_PATTERN, re.DOTALL)

SVG_DIM_PATTERN = r"width='(\d+px)'\s+height='(\d+px)"
SVG_DIM_REGEX = re.compile(SVG_DIM_PATTERN)

HEX_PATTERN = r'^#([A-Fa-f0-9]{6}|[A-Fa-f0-9]{3})$'
HEX_REGEX = re.compile(HEX_PATTERN)


def naive_svg_check(svg_string):
    """Validate SVG format (naive)."""
    return SVG_REGEX.match(svg_string) is not None


def svg_dimensions(svg_string):
    """Return dimensions of an (rdkit) SVG string."""
    matches = SVG_DIM_REGEX.findall(svg_string)
    if not matches:
        return (None, None)
    dims = map(lambda x: int(x.replace('px', '')), matches[0])
    return tuple(dims)


def insert_random_node_attribute(graph, key, high=1, low=0):
    """Add a random attribute to nodes in a graph."""
    for _, data in graph.nodes(data=True):
        value = random.uniform(low, high)
        data[key] = value


def is_valid_hex(hex):
    """Validate hexadecimal color code."""
    if hex is None:
        return False
    if HEX_REGEX.search(hex):
        return True
    return False


def test_smiles_to_svg():
    smi = 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'
    img = vis_utils.smiles_to_svg(smi)  # smiles to SVG
    assert img is not None
    assert naive_svg_check(img) is True
    # Check size updates.
    dims = (450, 400)
    img = vis_utils.smiles_to_svg(smi, size=dims)
    assert svg_dimensions(img) == dims
    # Check drawing options (clear background).
    drawOpts = rdMolDraw2D.MolDrawOptions()
    drawOpts.clearBackground = True
    img = vis_utils.smiles_to_svg(smi, draw_options=drawOpts)
    assert '</rect>' in img  # <rect> exists with background
    drawOpts.clearBackground = False
    img = vis_utils.smiles_to_svg(smi, draw_options=drawOpts)
    assert '</rect>' not in img


@suppress_rdlogger()
def test_smiles_to_image():
    # These aren't paticularly great tests for this function...
    smi = 'Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1'
    img = vis_utils.smiles_to_image(smi)  # smiles to SVG
    assert img is not None
    assert img != 'data:image/svg+xml;charset=utf-8,'
    null_smi = 'xxx'
    img = vis_utils.smiles_to_image(null_smi)
    assert img == 'data:image/svg+xml;charset=utf-8,'


def test_embed_node_mol_images(network):
    # Embed images into node attributes.
    vis_utils.embed_node_mol_images(network)
    for _, data in network.nodes(data=True):
        img = data.get('img', None)
        assert img is not None
    # Remove images from node attributes.
    vis_utils.remove_node_mol_images(network)
    for _, data in network.nodes(data=True):
        img = data.get('img', None)
        assert img is None


def test_color_nodes_by_attribute(network):
    key = 'attr'
    insert_random_node_attribute(network, key)
    # Color scaffold nodes.
    vis_utils.color_scaffold_nodes_by_attribute(network, key, 'BuPu')
    for _, data in network.get_scaffold_nodes(data=True):
        c = data.get('color', None)
        assert c is not None
        assert is_valid_hex(c)
    # Color molecule nodes.
    cmap = plt.get_cmap('hot')
    vis_utils.color_molecule_nodes_by_attribute(network, key, cmap, 'col')
    for _, data in network.get_molecule_nodes(data=True):
        c = data.get('col', None)
        assert c is not None
        assert is_valid_hex(c)


def test_root_node(network):
    vis_utils.add_root_node(network)
    assert network.has_node('root') is True
    assert network.in_degree('root') == 0
    vis_utils.remove_root_node(network)
    assert network.has_node('root') is False
