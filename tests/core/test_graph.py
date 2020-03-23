"""
scaffoldgraph tests.core.test_graph
"""

from rdkit import Chem

from scaffoldgraph.core.graph import *


def test_init_molecule_name():
    x = Chem.MolFromSmiles('CCC')
    assert bool(x.HasProp('_Name')) is False
    init_molecule_name(x)
    assert x.HasProp('_Name')
    assert x.GetProp('_Name') is not None
    assert x.GetProp('_Name') != ''


def test_graph_subclass():
    assert issubclass(ScaffoldGraph, nx.DiGraph)
    assert issubclass(ScaffoldGraph, ABC)
