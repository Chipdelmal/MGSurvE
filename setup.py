#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import setuptools
from version import version as this_version


this_directory =  os.path.abspath(os.path.dirname(__file__))
version_path = os.path.join(this_directory, 'MGSurvE', '_version.py')
print(this_version)
with open(version_path, 'wt') as fversion:
    fversion.write('__version__ = "'+this_version+'"')

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name='MGSurvE',
    version=this_version,
    url='https://github.com/Chipdelmal/MGSurvE',
    author='Hector M. Sanchez C.',
    author_email='sanchez.hmsc@berkeley.edu',
    description="MGSurvE",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    install_requires=[
        'deap', 'numpy', 'scikit-learn', 'scipy', 'matplotlib==3.5.2', 'sympy',
        'ipython', 'jupyter', 'pandas', 'compress-pickle', 'dill', 
        'vincenty', 'haversine', 'networkx', 'pointpats', 'libpysal',
        'geopandas'
    ],
    extras_require={
        'dev': [
            'pytest', 'ipykernel'
            'twine', 'wheel', 
            'sphinx', 'sphinx_rtd_theme'
        ],
    },
    license='GPLv3',
    classifiers=[
        'Programming Language :: Python :: 3.9',
        "Programming Language :: Python :: 3.10",
        "Operating System :: OS Independent",
    ]
 )