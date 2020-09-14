"""
scaffoldgraph.io

Contains functions for reading molecules from various formats.
    - SMILES
    - SDF
    - DataFrame

Constains function for writing ScaffoldGraphs to various formats.
    - SDF
    - TSV
"""

from .dataframe import read_dataframe
from .sdf import read_sdf
from .smiles import read_smiles_file

__all__ = ['read_sdf', 'read_smiles_file', 'read_dataframe']
