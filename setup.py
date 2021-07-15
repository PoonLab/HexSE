#!/usr/bin/env python3
from distutils.core import setup
from setuptools.command.install import install


setup(
    name='HexSE',
    version="0.0.1",
    description='Simulating evolution in overlapping reading frames',
    packages=['src'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'scipy',
        'numpy',
        'biopython',
        'yaml',
        'datetime',
        'random',
        'copy',
        'sys'
    ],
    
    python_requires='>=3.6'
)
