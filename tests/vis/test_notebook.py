"""
scaffoldgraph tests.vis.test_notebook
"""

from pytest import mark


@mark.filterwarnings('ignore::UserWarning')
def test_resources():
    # Import inside function to supress user warning.
    from scaffoldgraph.vis.notebook import cytoscape
    check = cytoscape.DEFAULT_STYLE
    assert check.parent.exists()  # resource directory
    assert check.exists()  # cytoscape.json
    style = cytoscape.read_style_file(str(check))
    assert style is not None
