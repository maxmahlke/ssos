#!/usr/bin/env python3
from setuptools import setup, find_packages
from os import path

setup(
    name='ssos',
    version='0.9.dev',
    description='ssos Pipeline',
    url='https://github.com/maxmahlke/ssos',
    author='Max Mahlke',
    author_email='max.mahlke@cab.inta-csic.es',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    entry_points={
            'console_scripts': ['ssos = ssos.__main__:main']
            }
    )