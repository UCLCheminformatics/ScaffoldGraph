"""
scaffoldgraph setup.py
"""

from setuptools import setup, find_packages
from pathlib import Path

__version__ = '1.0.4'

url = 'https://github.com/UCLCheminformatics/scaffoldgraph'

description = 'ScaffoldGraph is an open-source cheminformatics library, built using RDKit and \
NetworkX for generating scaffold networks and scaffold trees.'

root = Path(__file__).parent.resolve()

requires_path = root / 'requirements.txt'
with requires_path.open('r', encoding='utf8') as f:
    install_requires = [line.strip() for line in f]
    install_requires.remove('rdkit')

readme_path = root / 'README.md'
with readme_path.open('r', encoding='utf-8') as f:
    long_description = f.read()

setup_requires = ['pytest-runner']
tests_require = ['pytest', 'pytest-cov']

entry_points = {
    'console_scripts': [
        'scaffoldgraph = scaffoldgraph.scripts.run:scaffoldgraph_main',
    ]
}

setup(
    name='ScaffoldGraph',
    version=__version__,
    description=description,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Oliver Scott',
    author_email='oliver.scott.17@ucl.ac.uk',
    url=url,
    download_url='{}/archive/{}.tar.gz'.format(url, __version__),
    license='MIT',
    keywords=[
        'rdkit',
        'networkx',
        'cheminformatics',
        'scaffolds',
        'scaffold tree',
        'scaffold network'
    ],
    install_requires=install_requires,
    setup_requires=setup_requires,
    tests_require=tests_require,
    entry_points=entry_points,
    packages=find_packages(
        exclude=['tests.*', 'tests']
    ),
)
