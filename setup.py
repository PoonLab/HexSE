#!/usr/bin/env python3
from distutils.core import setup
from setuptools.command.install import install


setup(
    name='ovrf',
    description='Simulation of molecular evolution with overlapping reading frames',
    packages=['src'],
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent'
    ],
    install_requires=[
        'scipy',
        'numpy',
        'biopython'
    ],
    python_requires='>=3.6'
)
