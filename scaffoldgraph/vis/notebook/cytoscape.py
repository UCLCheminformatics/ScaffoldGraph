"""
scaffoldgraph.vis.notebook.cytoscape
"""

import warnings
import json

from pathlib import Path

from scaffoldgraph.vis.base import Visualizer
from scaffoldgraph.vis.utils import embed_node_mol_images

try:
    import ipycytoscape as cy
    _cytoscape_available = True
except ImportError:
    _cytoscape_available = False
    warnings.warn(
        'ipycytoscape could not be imported and is required '
        'for generating cytoscape based visualizations.'
    )

DEFAULT_STYLE = Path(__file__).parent.resolve() / 'resources' / 'cytoscape.json'

DEFAULT_LAYOUT = {
    'name': 'dagre',
    'nodeSpacing': 50,
    'edgeLengthVal': 50
}


def read_style_file(path):
    """Read a JSON style file (cytoscape).

    Parameters
    ----------
    path : str
        File path to style file.

    Returns
    -------
    style : dict
        Style dictionary.

    """
    with open(path, 'r') as f:
        style = json.load(f)
    return style


class CytoscapeVisualizer(Visualizer):
    """Class for creating visualizations using ipycytoscape.

    This visualizer renders scaffoldgraphs as interactive
    networks using cytoscape. The visualizer is flexible
    allowing users to customize the output defining the
    style and layout options.

    Notes
    -----
    visualizer is intended to be used within a jupyter notebook.

    ipycytoscape must be installed to use this feature.

    The code for this feature was inspired/adpated from:
    .. _Blogpost: https://iwatobipen.wordpress.com/2020/03/30/draw-scaffold-tree
    -as-network-with-molecular-image-rdkit-cytoscape/

    Examples
    --------
    Create a visualization for a whole graph.

    >>> from scaffoldgraph.vis.notebook import cytoscape
    >>> import scaffoldgraph as sg
    >>> tree = sg.ScaffoldTree.from_sdf('my_sdf.sdf')
    >>> visualizer = cytoscape.CytoscapeVisualizer(tree)
    >>> visualizer.draw()

    Use a different layout.

    >>> visualizer.draw(layout_kwargs={'name': 'breadthfirst'})

    Draw a subgraph starting from a molecule node.

    >>> visualizer.draw_for_molecule('CHEMBL1997663')

    Draw a subgraph starting from a scaffold node.

    >>> visualizer.draw_for_scaffold('c1ccc(CNc2ccccc2)cc1')

    """
    def __init__(
            self,
            graph,
            style=None,
            refresh_images=False,
            rd_draw_options=None,
            mol_img_size=(350, 300),
    ):
        """Initialize the cytoscape visualizer.

        Parameters
        ----------
        graph : ScaffoldGraph
            A ScaffoldGraph object to draw.
        style : list, optional
            A list of dicts specifying the style to pass
            to the cytoscape widget, for more details
            see the ipycytoscape documentation. If None
            a default style is used and can be updated
            after initialization.
        refresh_images: bool, optional
            If True remove all embeded images from the
            input graph and regenerate when required.
            The default is False.
        rd_draw_options: rdkit.Chem.Draw.rdMolDraw2D.MolDrawOptions, optional
            Specify options for molecule drawing. Requires a
            `MolDrawOptions` object or `None`.
            The default is None.
        mol_img_size: tuple, optional
            Specify the size of the node images. Format is
            `(width, height)`. Note that if changed from
            default the style will have to be updated.
            The default is `(350, 300)`.

        """
        super(CytoscapeVisualizer, self).__init__(
            graph,
            requires_tree=False,
            refresh_images=refresh_images,
        )
        self._drawopts = rd_draw_options
        self._img_size = mol_img_size
        self._style = style if style else read_style_file(DEFAULT_STYLE)

    @property
    def style(self):
        """list : returns the cytoscape style associated."""
        return self._style

    @style.setter
    def style(self, style):
        assert isinstance(style, list),\
            'style must be a list object'
        self._style = style

    @staticmethod
    def _cytoscape_validate():
        if _cytoscape_available is False:
            raise RuntimeError('ipycytoscape is not available')

    def _draw(self, subgraph, layout_kwargs):
        """Private: create the cytoscape widget from a subgraph."""
        if subgraph.number_of_nodes() >= 100:
            warnings.warn('graphs with > 100 nodes may be slow to render')
        embed_node_mol_images(
            subgraph,
            size=self._img_size,
            draw_options=self._drawopts,
        )
        layout = {}
        layout.update(DEFAULT_LAYOUT)
        if layout_kwargs:
            layout.update(layout_kwargs)
        widget = cy.CytoscapeWidget()
        widget.set_style(self._style)
        widget.set_layout(**layout)
        widget.graph.add_graph_from_networkx(
            subgraph, directed=True
        )
        return widget

    def draw(self, layout_kwargs=None):
        """Draw the entire scaffoldgraph.

        Parameters
        ----------
        layout_kwargs : dict, optional
            arguments to pass to the CytoscapeWidget.set_layout
            function.

        Returns
        -------
        widget : ipycytoscape.CytoscapeWidget

        """
        self._cytoscape_validate()
        return self._draw(self._graph, layout_kwargs)

    def draw_for_molecule(self, molecule_id, layout_kwargs=None):
        """Draw subgraph starting from a query molecule.

        Parameters
        ----------
        molecule_id : str
            Molecule node identifier.
        layout_kwargs : dict, optional
            arguments to pass to the CytoscapeWidget.set_layout
            function.

        Returns
        -------
        widget : ipycytoscape.CytoscapeWidget

        """
        self._cytoscape_validate()
        subgraph = self._subgraph_from_mol(molecule_id)
        return self._draw(subgraph, layout_kwargs)

    def draw_for_scaffold(self, scaffold_id, traversal='child', layout_kwargs=None):
        """Draw subgraph starting from a query scaffold.

        Parameters
        ----------
        scaffold_id : str
            Scaffold node identifier.
        traversal : str {'parent', 'child', 'bidirectional'}
            The direction of traversal to create the subgraph.
            If 'bidirectional' both directions are considered.
            The default is 'child'.
        layout_kwargs : dict, optional
            arguments to pass to the CytoscapeWidget.set_layout
            function.

        Returns
        -------
        widget : ipycytoscape.CytoscapeWidget

        """
        self._cytoscape_validate()
        subgraph = self._subgraph_from_scf(scaffold_id, traversal)
        return self._draw(subgraph, layout_kwargs)
