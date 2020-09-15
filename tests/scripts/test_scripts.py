"""
scaffoldgraph tests.scripts.test_scripts
"""

import tempfile
import pathlib
import pytest
import os

from subprocess import Popen, PIPE


TEST_DATA_DIR = pathlib.Path(__file__).resolve().parent / '..' / 'data'


def check_generate_structure(fn):
    with open(fn, 'r') as f:
        lines = f.readlines()
    headings = lines[0].strip().split('\t')
    assert 'HIERARCHY' in headings
    assert 'SMILES' in headings
    assert 'SUBSCAFFOLDS' in headings
    assert 'MOLECULES' in headings
    assert 'ANNOTATIONS' in headings
    assert len(lines) > 1


def check_aggregate_structure(fn):
    with open(fn, 'r') as f:
        lines = f.readlines()
    headings = lines[0].strip().split('\t')
    assert 'HIERARCHY' in headings
    assert 'SMILES' in headings
    assert 'SUBSCAFFOLDS' in headings
    assert 'ID' in headings
    assert len(lines) > 1
    smiles = lines[-1].strip().split('\t')[2]
    return smiles


def check_select_structure(fn):
    with open(fn, 'r') as f:
        lines = f.readlines()
    headings = lines[0].strip().split('\t')
    assert 'HIERARCHY' in headings
    assert 'SMILES' in headings
    assert 'SUBSCAFFOLDS' in headings
    assert 'ID' in headings
    assert len(lines) > 1


# test all utilities in one
# skip: pytest -m "not slow"
@pytest.mark.slow
def test_cli():
    funcs = ['tree', 'network', 'hiers']
    fn = str(TEST_DATA_DIR / 'test_smiles.smi')
    with tempfile.TemporaryDirectory() as tmp:
        for func in funcs:

            # Test graph generation
            out = os.path.join(tmp, 'output.tmp')
            args = ['scaffoldgraph', func, fn, out]
            p2 = Popen(args, stdout=PIPE, stderr=PIPE)
            stdout, _ = p2.communicate()
            assert stdout is not None
            assert os.path.exists(out)
            check_generate_structure(out)

            # Test graph aggregation
            out2 = os.path.join(tmp, 'output.txt')
            args = [
                'scaffoldgraph', 'aggregate', out, out2
            ]
            p2 = Popen(args, stdout=PIPE, stderr=PIPE)
            stdout, _ = p2.communicate()
            assert stdout is not None
            assert os.path.exists(out)
            smiles = check_aggregate_structure(out2)

            # Test graph selection
            test_smi = os.path.join(tmp, 'test.smi')
            out3 = os.path.join(tmp, 'select.txt')
            with open(test_smi, 'w') as smi:
                smi.write(f'{smiles} fake_scaffold_id')
            args = [
                'scaffoldgraph', 'select', out2, test_smi, out3
            ]
            p2 = Popen(args, stdout=PIPE, stderr=PIPE)
            stdout, _ = p2.communicate()
            assert stdout is not None
            assert os.path.exists(out3)
            check_select_structure(out3)
