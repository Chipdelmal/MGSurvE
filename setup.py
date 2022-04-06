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
        'deap', 'numpy', 'matplotlib', 'ipython',
        'jupyter', 'pandas', 'compress-pickle',
        'scikit-learn', 'scipy', 'vincenty',
        'pytest', 'networkx', 'pointpats', 'geopandas',
        'haversine', 'dill', 'pytest'
    ],
    license='MIT',
    classifiers=[
        "Programming Language :: Python :: 3.9",
        "Operating System :: OS Independent",
    ]
 )