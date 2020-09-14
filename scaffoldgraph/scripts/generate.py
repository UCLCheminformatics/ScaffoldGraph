"""
scaffoldgraph.scripts.generate
"""

import datetime
import time

from loguru import logger

from scaffoldgraph import ScaffoldNetwork, ScaffoldTree, HierS
from scaffoldgraph.prioritization import ScaffoldRuleSet
from scaffoldgraph.io import tsv

from .misc import file_format

start_message = """
Running ScaffoldGraph ({command}) Generation with options:
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


def _get_graph_cls(name):
    """Get scaffoldgraph class from name string."""
    if name == 'network':
        return ScaffoldNetwork
    elif name == 'tree':
        return ScaffoldTree
    elif name == 'hiers':
        return HierS
    else:
        msg = f'scaffold graph type: {name} not known'
        raise ValueError(msg)


def _maybe_ruleset(args):
    """Return a ScaffoldRuleset if specified in CLI arguments."""
    ruleset = None
    if 'ruleset' in args and args.ruleset is not None:
        filename = args.ruleset
        ruleset = ScaffoldRuleSet.from_rule_file(filename)
    return ruleset


def generate_cli(args):
    """Run scaffoldgraph generation for CLI utility."""
    graph_cls = _get_graph_cls(args.command)
    graph_name = graph_cls.__name__
    ruleset = _maybe_ruleset(args)

    if not args.silent:
        print(
            start_message.format(
                command=graph_name,
                input=args.input,
                output=args.output,
                max_r=args.max_rings
            )
        )

    logger.info(f'Generating {graph_name} Graph...')
    fmt, zipped = file_format(args.input)
    start = time.time()

    if fmt == 'SDF':
        sg = graph_cls.from_sdf(
            args.input,
            ring_cutoff=args.max_rings,
            progress=args.silent is False,
            zipped=zipped,
            prioritization_rules=ruleset,
        )
    elif fmt == 'SMI':
        sg = graph_cls.from_smiles_file(
            args.input,
            ring_cutoff=args.max_rings,
            progress=args.silent is False,
            prioritization_rules=ruleset,
        )
    else:
        raise ValueError('input file format is not currently supported')

    tsv.write_tsv(sg, args.output, write_ids=False)
    logger.info(f'{graph_name} Graph Generation Complete...')
    elapsed = datetime.timedelta(seconds=round(time.time() - start))

    if not args.silent:
        print(
            stop_message.format(
                molecules=sg.num_molecule_nodes,
                scaffolds=sg.num_scaffold_nodes,
                time=elapsed,
                output=args.output
            )
        )
