#!/usr/bin/env python3
from setuptools import setup, find_packages

with open("README.md", "r") as rm:
    long_description = rm.read()

setup(
    name="ssos",
    version="1.3.5",
    description=f"The ssos Pipeline - Detection of Solar System Objects "
    f"in astronomical images",
    url="https://github.com/maxmahlke/ssos",
    author="Max Mahlke",
    author_email="max.mahlke@oca.eu",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=find_packages(exclude=["contrib", "docs", "tests"]),
    keywords=["astronomy astrophysics solar-system data pipeline"],
    setup_requires=[
        "astropy",
        "pandas>=0.23.0",
        "matplotlib",
        "numpy",
        "scipy",
        "statsmodels>=0.9.0",
        "tqdm>=4.40.2",
    ],
    install_requires=[
        "astropy",
        "pandas>=0.23.0",
        "matplotlib",
        "numpy",
        "scipy",
        "statsmodels>=0.9.0",
        "tqdm>=4.40.2",
    ],
    entry_points={"console_scripts": ["ssos = ssos.__main__:main"]},
    include_package_data=True,
    python_requires=">=3.6",
)
