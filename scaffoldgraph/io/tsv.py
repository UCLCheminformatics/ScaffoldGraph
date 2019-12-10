"""
scaffoldgraph.io.tsv
"""

import csv


def write_tsv(scaffold_graph, output_file, write_ids=False):
    """Write a graph to a file in a TSV format.

    Used by the CLI utility.

    Parameters
    ----------
    scaffold_graph : sg.ScaffoldGraph
        An scaffold graph to write to a file
    output_file : str
        Path to output file
    """

    N = scaffold_graph.num_scaffold_nodes
    sorted_scaffolds = sorted(scaffold_graph.get_scaffold_nodes(data=True), key=lambda x: x[1]['hierarchy'])

    if write_ids:
        field_names = ['ID', 'HIERARCHY', 'SMILES', 'SUBSCAFFOLDS']
        mapping = dict(zip([s[0] for s in sorted_scaffolds], range(0, N)))
    else:
        field_names = ['HIERARCHY', 'SMILES', 'SUBSCAFFOLDS', 'MOLECULES', 'ANNOTATIONS']
        mapping = None

    with open(output_file, 'w') as output:

        writer = csv.DictWriter(output, delimiter='\t', fieldnames=field_names)
        writer.writeheader()

        for node, data in sorted_scaffolds:
            line = dict.fromkeys(field_names)
            line['SMILES'] = node
            line['HIERARCHY'] = data['hierarchy']

            subscaffolds = list(scaffold_graph.predecessors(node))
            if write_ids:
                line['SUBSCAFFOLDS'] = ', '.join([str(mapping[s]) for s in subscaffolds])
                line['ID'] = str(mapping[node])
            else:
                line['SUBSCAFFOLDS'] = ', '.join(subscaffolds)
                ancestors = scaffold_graph.successors(node)
                molecules, annotations = [], set()
                for a in ancestors:
                    try:
                        if scaffold_graph.nodes[a]['type'] == 'molecule':
                            molecules.append(a)
                            edge = scaffold_graph.edges[(node, a)]
                            annotations.add(edge['annotation'])
                    except KeyError:
                        continue
                line['MOLECULES'] = ', '.join(molecules)
                line['ANNOTATIONS'] = ', '.join(annotations)

            writer.writerow(line)
