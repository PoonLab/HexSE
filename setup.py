#!/usr/bin/env python3
from distutils.core import setup
from setuptools.command.install import install


setup(
    name='HexSE',
    version="0.0.1",
    description='Simulating evolution in overlapping reading frames',
    packages=['hexse'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'scipy',
        'numpy',
        'biopython',
        'pyyaml', 
        'tqdm'
    ],
    
    python_requires='>=3.6'
)
