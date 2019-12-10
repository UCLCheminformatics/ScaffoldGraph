"""
scaffoldgraph.scripts.operations
"""

import datetime
import time

from loguru import logger
from rdkit import RDLogger
from rdkit.Chem import SDWriter, MolFromSmiles, MolToSmiles

from .misc import file_format
from ..core import get_murcko_scaffold
from ..io import smiles, sdf

rdlogger = RDLogger.logger()


class ScaffoldRecord(object):
    """Class holds a scaffold read from a TSV output file"""
    __slots__ = (
        'id',
        'hierarchy',
        'smiles',
        'subscaffolds',
        'molecules',
        'annotations'
    )

    def __init__(self):
        self.id = None
        self.hierarchy = 0
        self.smiles = None
        self.subscaffolds = []
        self.molecules = set()
        self.annotations = set()

    def add_subscaffold(self, record):
        self.subscaffolds.append(record)

    def add_molecule(self, molecule):
        self.molecules.add(molecule)

    def add_annotation(self, annotation):
        self.annotations.add(annotation)


class ScaffoldFileIterator(object):
    """Read an output TSV file into ScaffoldRecords"""

    def __init__(self, fw, reverse=False):
        self._fw = fw
        header = fw.readline()
        self.columns = header.strip().split('\t')
        self.column_map = self.map_columns(self.columns)
        if reverse is True:
            self._fw = reversed(self._fw.readlines())

    def map_columns(self, columns):
        return {col.upper(): idx for (idx, col) in enumerate(columns)}

    def process_record(self, record):
        s = ScaffoldRecord()  # empty record

        idx = self.column_map.get('ID', None)
        s.id = int(record[idx]) if idx is not None else None

        idx = self.column_map.get('HIERARCHY', None)
        s.hierarchy = int(record[idx]) if idx is not None else 0

        idx = self.column_map.get('SMILES', None)
        s.smiles = str(record[idx]).strip() if idx is not None else None

        idx = self.column_map.get('SUBSCAFFOLDS', None)
        val = record[idx].strip().split(',') if idx is not None else []
        for subscaffold in val:
            if subscaffold != '':
                sub = ScaffoldRecord()
                try:
                    sub.id = int(subscaffold)
                except ValueError:
                    sub.smiles = str(subscaffold).strip()
                s.add_subscaffold(sub)

        idx = self.column_map.get('MOLECULES', None)
        val = record[idx].strip().split(',') if idx is not None else []
        for molecule in val:
            if molecule != '':
                s.add_molecule(molecule.strip())

        idx = self.column_map.get('ANNOTATIONS', None)
        val = record[idx].strip().split(',') if idx is not None else []
        for annotation in val:
            if annotation != '':
                s.add_annotation(annotation.strip())

        return s

    def __iter__(self):
        return self

    def __next__(self):
        line = next(self._fw)
        record = line.strip('\n').split('\t')
        return self.process_record(record)


class AggregateCLI(object):
    """Aggregate output TSV files (CLI)"""

    def __init__(self, args):

        self.args = args
        self.inputs = args.input

        if args.sdf:
            rdlogger.setLevel(4)
            self.output = SDWriter(args.output)
        else:
            self.output = open(args.output, 'w')

        self.mol_map = open(args.map_mols, 'w') if args.map_mols else None
        if self.mol_map:
            self.mol_map.write('MOLECULE_ID\tSCAFFOLD_ID\n')

        self.ann_map = open(args.map_annotations, 'w') if args.map_annotations else None
        if self.ann_map:
            self.ann_map.write('SCAFFOLD_ID\tANNOTATIONS\n')

        self.current_id = 0
        self.duplicates = 0
        self.table = {}

    def aggregate(self):
        if not self.args.sdf:
            self.output.write('ID\tHIERARCHY\tSMILES\tSUBSCAFFOLDS\n')
        for file in self.inputs:
            logger.info(f'Processing file: {file}...')
            with open(file, 'r') as fw:
                self.process_file(fw)

    def process_file(self, file):
        reader = ScaffoldFileIterator(file)
        for scaffold in reader:
            s_smiles = scaffold.smiles
            write = False
            if s_smiles in self.table:
                scaffold.id = self.table[s_smiles]['ID']
            else:
                scaffold.id = self.current_id
                self.table[s_smiles] = dict(ID=self.current_id, PARENTS=[])
                self.current_id += 1
                write = True
            missing = []
            for idx, parent in enumerate(scaffold.subscaffolds):
                p_smiles = parent.smiles
                if p_smiles in self.table:
                    parent.id = self.table[p_smiles]['ID']
                else:
                    missing.append(idx)
            for m in sorted(missing, reverse=True):
                del scaffold.subscaffolds[m]
            if write is True:
                self.write_scaffold(scaffold)
                self.write_extra_outputs(scaffold)
            else:
                self.duplicates += 1
                self.write_extra_outputs(scaffold)

    def write_scaffold(self, scaffold):
        subscaffolds = ', '.join([str(s.id) for s in scaffold.subscaffolds])
        if self.args.sdf:
            molecule = MolFromSmiles(scaffold.smiles)
            if molecule is not None:
                molecule.SetProp('_Name', str(scaffold.id))
                molecule.SetIntProp('HIERARCHY', scaffold.hierarchy)
                molecule.SetProp('SMILES', scaffold.smiles)
                molecule.SetProp('SUBSCAFFOLDS', subscaffolds)
                self.output.write(molecule)
            else:
                logger.warning(f'Failed to parse scaffold: {scaffold.smiles}')
        else:
            self.output.write('{0}\t{1}\t{2}\t{3}\n'.format(
                scaffold.id,
                scaffold.hierarchy,
                scaffold.smiles,
                subscaffolds))

    def write_extra_outputs(self, scaffold):
        # Write molecule --> scaffold ID file
        if self.mol_map is not None:
            for molecule in scaffold.molecules:
                self.mol_map.write('{0}\t{1}\n'.format(
                    molecule, scaffold.id))
        # Write scaffold ID --> annotation file
        if self.ann_map is not None:
            for annotation in scaffold.annotations:
                self.ann_map.write('{0}\t{1}\n'.format(scaffold.id, annotation))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.output.close()
        if self.mol_map is not None:
            self.mol_map.close()
        if self.ann_map is not None:
            self.ann_map.close()


class SelectCLI(object):
    """Select scaffolds using a molecular query on an aggregated TSV (CLI)."""

    def __init__(self, args):

        self.args = args
        self.q_input = args.input_query
        self.g_input = open(args.input_graph, 'r')

        if args.sdf:
            rdlogger.setLevel(4)
            self.output = SDWriter(args.output)
        else:
            self.output = open(args.output, 'w')

        self.query = set()
        self.matching_parents = set()
        self.count = 0

    def select(self):
        if not self.args.sdf:
            self.output.write('ID\tHIERARCHY\tSMILES\tSUBSCAFFOLDS\n')
        self.load_query()
        logger.info('Processing query...')
        reader = ScaffoldFileIterator(self.g_input, reverse=True)
        for scaffold in reader:
            match = False
            if scaffold.smiles in self.query:
                match = True
            if scaffold.id in self.matching_parents:
                match = True
            if match is True:
                self.count += 1
                self.write_scaffold(scaffold)
                for s in scaffold.subscaffolds:
                    self.matching_parents.add(s.id)

    def load_query(self):
        logger.info('Reading molecular query...')
        file = None
        fmt = file_format(self.q_input)
        if fmt == 'SMI':
            supplier = smiles.read_smiles_file(self.q_input)
        elif fmt == 'SDF':
            rdlogger.setLevel(4)
            file = open(self.q_input, 'rb')
            supplier = sdf.read_sdf(file)
        else:
            raise ValueError('input file format not currently supported')
        for molecule in supplier:
            if molecule is not None:
                s = get_murcko_scaffold(molecule)
                self.query.add(MolToSmiles(s))
        if file is not None:
            file.close()
        logger.info(f'Read {len(self.query)} query scaffolds')

    def write_scaffold(self, scaffold):
        subscaffolds = ', '.join([str(s.id) for s in scaffold.subscaffolds])
        if self.args.sdf:
            molecule = MolFromSmiles(scaffold.smiles)
            if molecule is not None:
                molecule.SetProp('_Name', str(scaffold.id))
                molecule.SetIntProp('HIERARCHY', scaffold.hierarchy)
                molecule.SetProp('SMILES', scaffold.smiles)
                molecule.SetProp('SUBSCAFFOLDS', subscaffolds)
                self.output.write(molecule)
            else:
                logger.warning(f'Failed to parse scaffold: {scaffold.smiles}')
        else:
            self.output.write('{0}\t{1}\t{2}\t{3}\n'.format(
                scaffold.id,
                scaffold.hierarchy,
                scaffold.smiles,
                subscaffolds))

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.g_input.close()
        self.output.close()


start_message = """
Running ScaffoldGraph {command} Operation with options:
    Input file(s):  {input}
    Output file:    {output}
"""

stop_message_agg = """
ScaffoldGraph {command} Complete:
    Scaffolds written: {scaffolds}
    Duplicates: {duplicates}
    Time elapsed: {time}

Output saved @ {output}
"""


def aggregate_cli(args):
    """Command line function for aggregating intermediate outputs"""

    if not args.silent:
        print(start_message.format(command='Aggregate', input=args.input,
                                   output=args.output))
    start = time.time()
    logger.info('Aggregating graphs...')
    with AggregateCLI(args) as aggregator:
        aggregator.aggregate()

    logger.info('Scaffold Graph Aggregation Complete.')
    elapsed = datetime.timedelta(seconds=round(time.time() - start))

    if not args.silent:
        print(stop_message_agg.format(command='Aggregate',
                                      scaffolds=aggregator.current_id,
                                      duplicates=aggregator.duplicates,
                                      time=elapsed, output=args.output))


stop_message_sel = """
ScaffoldGraph {command} Complete:
    Scaffolds written: {scaffolds}
    Time elapsed: {time}

Output saved @ {output}
"""


def select_cli(args):
    """Command line function for selecting a subset using a molecular query"""

    if not args.silent:
        print(start_message.format(command='Select',
                                   input=[args.input_graph, args.input_query],
                                   output=args.output))
    start = time.time()
    with SelectCLI(args) as selector:
        selector.select()

    logger.info('Scaffold Graph Selection Operation Complete.')
    elapsed = datetime.timedelta(seconds=round(time.time() - start))

    if not args.silent:
        print(stop_message_sel.format(command='Select',
                                      scaffolds=selector.count,
                                      time=elapsed, output=args.output))
