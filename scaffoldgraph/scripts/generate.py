"""
scaffoldgraph.scripts.generate
"""

import datetime
import time

from loguru import logger

from .misc import file_format
from .. import ScaffoldNetwork, ScaffoldTree, HierS
from ..io import tsv

start_message = """
Running ScaffoldGraph {command} Generation with options:
    Input file:     {input}
    Output file:    {output}
    Maximum rings:  {max_r}
"""

stop_message = """
ScaffoldGraph Generation Complete:
    Molecules written: {molecules}
    Scaffolds written: {scaffolds}
    Time elapsed: {time}
    
Output saved @ {output}
"""


def network_cli(args):
    """Run network generation for CLI utility"""
    if not args.silent:
        print(start_message.format(command='Network', input=args.input,
                                   output=args.output, max_r=args.max_rings))
    logger.info('Generating Scaffold Network...')
    fmt = file_format(args.input)
    start = time.time()
    if fmt == 'SDF':
        sg = ScaffoldNetwork.from_sdf(args.input, args.max_rings, progress=args.silent is False)
    elif fmt == 'SMI':
        sg = ScaffoldNetwork.from_smiles_file(args.input, args.max_rings, progress=args.silent is False)
    else:
        raise ValueError('input file format not currently supported')

    tsv.write_tsv(sg, args.output, write_ids=False)

    logger.info('Scaffold Network Generation Complete...')
    elapsed = datetime.timedelta(seconds=round(time.time() - start))
    if not args.silent:
        print(stop_message.format(molecules=sg.num_molecule_nodes, scaffolds=sg.num_scaffold_nodes,
                                  time=elapsed, output=args.output))


def hiers_cli(args):
    """Run HierS network generation for CLI utility"""
    if not args.silent:
        print(start_message.format(command='HierS', input=args.input,
                                   output=args.output, max_r=args.max_rings))
    logger.info('Generating HierS Scaffold Network...')
    fmt = file_format(args.input)
    start = time.time()
    if fmt == 'SDF':
        sg = HierS.from_sdf(args.input, args.max_rings, progress=args.silent is False)
    elif fmt == 'SMI':
        sg = HierS.from_smiles_file(args.input, args.max_rings, progress=args.silent is False)
    else:
        raise ValueError('input file format not currently supported')

    tsv.write_tsv(sg, args.output, write_ids=False)

    logger.info('HierS Scaffold Network Generation Complete...')
    elapsed = datetime.timedelta(seconds=round(time.time() - start))
    if not args.silent:
        print(stop_message.format(molecules=sg.num_molecule_nodes, scaffolds=sg.num_scaffold_nodes,
                                  time=elapsed, output=args.output))


def tree_cli(args):
    """Run tree generation for CLI utility"""
    if not args.silent:
        print(start_message.format(command='Tree', input=args.input,
                                   output=args.output, max_r=args.max_rings))
    logger.info('Generating Scaffold Tree...')
    fmt = file_format(args.input)
    start = time.time()
    if fmt == 'SDF':
        sg = ScaffoldTree.from_sdf(args.input, args.max_rings, progress=args.silent is False)
    elif fmt == 'SMI':
        sg = ScaffoldTree.from_smiles_file(args.input, args.max_rings, progress=args.silent is False)
    else:
        raise ValueError('input file format not currently supported')

    tsv.write_tsv(sg, args.output, write_ids=False)

    logger.info('Scaffold Tree Generation Complete...')
    elapsed = datetime.timedelta(seconds=round(time.time() - start))
    if not args.silent:
        print(stop_message.format(molecules=sg.num_molecule_nodes, scaffolds=sg.num_scaffold_nodes,
                                  time=elapsed, output=args.output))
