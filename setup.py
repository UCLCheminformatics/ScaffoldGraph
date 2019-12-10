"""
scaffoldgraph setup
"""

from setuptools import setup, find_packages

__version__ = '0.1.0'
url = 'https://github.com/OliverBScott/scaffoldgraph'

description = 'ScaffoldGraph is an open-source cheminformatics library, built using RDKit and \
NetworkX for generating scaffold networks and scaffold trees.'

with open('requirements.txt') as f:
    install_requires = [l.strip() for l in f]
    install_requires.remove('rdkit')

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
    ],
    install_requires=install_requires,
    setup_requires=setup_requires,
    tests_require=tests_require,
    entry_points=entry_points,
    packages=find_packages(),
)
