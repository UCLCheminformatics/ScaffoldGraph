.. _installation:

Installation
============

To get started with ScaffoldGraph, install with pip:

.. code-block:: bash

    pip install scaffoldgraph

or with conda:

.. code-block:: bash

    conda config --add channels conda-forge
    conda install -c uclcheminformatics scaffoldgraph

or from source:

.. code-block:: bash

    wget https://github.com/UCLCheminformatics/ScaffoldGraph/archive/1.0.4.tar.gz
    tar -zxvf 1.0.4.tar.gz
    cd ScaffoldGraph-1.0.4
    python setup.py install

.. note:: When installing ScaffoldGraph with pip or from source, RDKit is assumed to be pre-installed as it cannot be installed using pip.
  For more details visit the RDKit `documentation <https://www.rdkit.org/docs/Install.html>`_.


Requirements
------------

ScaffoldGraph depends on the following packages:

- Python (3.5 +)
- RDKit_ (Cheminformatics toolkit)
- NetworkX_ (Graph analysis library)
- SciPy_ (Scientific computing library)
- tqdm_ (Progress bars)
- loguru_ (Logging library)

.. _RDKit: https://www.rdkit.org/docs/Overview.html/
.. _NetworkX: https://networkx.org/
.. _SciPy: https://www.scipy.org/
.. _tqdm: https://tqdm.github.io/
.. _loguru: https://loguru.readthedocs.io/en/stable/index.html


Optional Requirements
---------------------

ScaffoldGraph contains visualization functions depending on packages that are not installed by default.

**CytoscapeVisualizer**:
  The cytoscape visualizer depends on the ipycytoscape_ package. It can be installed with pip. For more details
  visit the ipycytoscape_ documentation.

.. code-block:: bash

    pip install ipycytoscape

or with conda:

.. code-block:: bash

    conda install -c conda-forge ipycytoscape


.. _ipycytoscape: https://ipycytoscape.readthedocs.io/en/latest/
