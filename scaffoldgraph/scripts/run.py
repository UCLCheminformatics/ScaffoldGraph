"""
scaffoldgraph.scripts.run

Module defines the CLI utility for ScaffoldGraph
"""

import argparse
import logging
import sys

from scaffoldgraph import __version__
from .generate import *
from .misc import TqdmHandler
from .operations import *

title = f"ScaffoldGraph {__version__}"
desc = "Generate Scaffold Networks and Scaffold Trees."

tqdm_format = "<light-green>scaffold-graph:</light-green> "
tqdm_format += "<yellow>{time:HH:mm:ss}</yellow> "
tqdm_format += "<red>{process}</red> "
tqdm_format += "<level>{level}:</level> "
tqdm_format += "{message}"

tqdm_handler = {
    'sink': TqdmHandler(logging.NOTSET),
    'format': tqdm_format,
    'level': 'INFO'
}

usage = 'scaffoldgraph <command> [<args>]'


def configure_logger(verbosity):
    """Configure the scaffoldgraph cli logger to use tqdm handler"""

    config = {'handlers': []}
    logger.enable('scaffoldgraph')

    if verbosity == 0:
        tqdm_handler['sink'].level = logging.CRITICAL
        tqdm_handler['level'] = 'CRITICAL'
    elif verbosity == 1:
        tqdm_handler['sink'].level = logging.ERROR
        tqdm_handler['level'] = 'ERROR'
    elif verbosity == 2:
        tqdm_handler['sink'].level = logging.WARNING
        tqdm_handler['level'] = 'WARNING'
    elif verbosity == 3:
        tqdm_handler['sink'].level = logging.INFO
        tqdm_handler['level'] = 'INFO'
    elif verbosity == 4:
        tqdm_handler['sink'].level = logging.DEBUG
        tqdm_handler['level'] = 'DEBUG'

    config["handlers"].append(tqdm_handler)
    logger.configure(**config)


def parent_parser():
    """Common arguments for all scaffoldgraph commands"""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('-v', '--verbosity', metavar='', type=int, default=3, choices=[0, 1, 2, 3, 4],
                        help='set logger verbosity [0, 1, 2, 3, 4] (default: 3)')
    parser.add_argument('-s', '--silent', action='store_true', help='silence console output (default: False)')
    return parser


def generate_parent_parser():
    """Creates a parent parser for generate commands (Network, Tree)."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument('input', help='input file (SDF, SMILES)')
    parser.add_argument('output', help='output file path')
    parser.add_argument('--max-rings', '-m', type=int, default=10, metavar='',
                        help='ignore molecules with # rings > (default: 10)')
    return parser


def scaffoldgraph_args():
    """ Defines CLI utility for ScaffoldGraph."""

    parser = argparse.ArgumentParser('scaffoldgraph', description=desc)
    parser.add_argument('--version', action='version', version=__version__)
    subparsers = parser.add_subparsers(title='command', dest='command')

    # network (generate a scaffold network from a SMILES or SDF file)
    network_parser = subparsers.add_parser('network', description='Generate a scaffold network',
                                           parents=[generate_parent_parser(), parent_parser()])
    network_parser.set_defaults(func=network_cli)

    # HierS (generate a HierS scaffold network from a SMILES or SDF file)
    hiers_parser = subparsers.add_parser('hiers', description='Generate a HierS type scaffold network',
                                         parents=[generate_parent_parser(), parent_parser()])
    hiers_parser.set_defaults(func=hiers_cli)

    # tree (generate a scaffold tree form a SMILES or SDF file)
    tree_parser = subparsers.add_parser('tree', description='Generate a scaffold tree',
                                        parents=[generate_parent_parser(), parent_parser()])
    tree_parser.set_defaults(func=tree_cli)

    # select (select a subgraph of a scaffold graph using a molecular query)
    select_parser = subparsers.add_parser('select', description="Select subgraph from a molecular query.",
                                          parents=[parent_parser()])
    select_parser.add_argument('input_graph', help='input aggregated graph file')
    select_parser.add_argument('input_query', help='input query file (SDF, SMILES)')
    select_parser.add_argument('output', help='output file path')
    select_parser.add_argument('-d', '--sdf', help='write output as an SDF', action='store_true')
    select_parser.set_defaults(func=select_cli)

    # aggregate (Aggregate intermediate scaffold graph files (TSV or PICKLE))
    aggregate_parser = subparsers.add_parser('aggregate', description='Aggregate scaffold graphs',
                                             parents=[parent_parser()])
    aggregate_parser.add_argument('input', nargs='+', help='input file(s) (TSV)')
    aggregate_parser.add_argument('output', help='output file path')
    aggregate_parser.add_argument('-m', '--map-mols', help='map molecule IDs from input to scaffold IDs, \
                                  and place result in given file', metavar='')
    aggregate_parser.add_argument('-a', '--map-annotations', help='map scaffold IDs to annotations, \
                                  and place result in given file', metavar='')
    aggregate_parser.add_argument('-d', '--sdf', help='write output as an SDF', action='store_true')
    aggregate_parser.set_defaults(func=aggregate_cli)

    return parser


def scaffoldgraph_main():
    """Run the CLI utility for ScaffoldGraph."""
    parser = scaffoldgraph_args()
    args = parser.parse_args(None if sys.argv[1:] else ['-h'])
    configure_logger(args.verbosity)
    try:
        args.func(args)
    except FileNotFoundError as e:
        logger.critical(f'Input file not found: {e.filename}')
    except ValueError as e:
        logger.critical(e)
    except RuntimeError as e:
        logger.critical(e)
    except MemoryError as e:
        logger.critical(e)
    except KeyboardInterrupt:
        logger.critical('scaffoldgraph process interrupted from keyboard')
    except Exception as e:
        logger.critical(f'Unknown error: {e}')
    finally:
        logger.info('Exiting scaffoldgraph...')
