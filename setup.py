#!/usr/bin/env python3
from setuptools import setup, find_packages
from os import path

setup(
    name='ssos',
    version='0.9.dev0',
    description='ssos Pipeline',
    url='https://github.com/maxmahlke/ssos',
    author='Max Mahlke',
    author_email='max.mahlke@cab.inta-csic.es',
    packages=find_packages(exclude=['contrib', 'docs', 'tests']),
    setup_requires=packages,
    install_requires = packages,
    entry_points={
            'console_scripts': ['ssos = ssos.__main__:main']
            }
    )