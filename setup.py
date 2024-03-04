#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Oct 1, 2021

@author: mwlkhoo
"""

from setuptools import setup, find_packages


def readme():
    with open("README.md") as f:
        return f.read()


setup(
    name="stats_utils",
    version="0.0.1",
    description="Statistics utilities for remoscope and corresponding paper",
    long_description=readme(),
    url="https://github.com/czbiohub-sf/remo-stats-utils",
    author="Bioengineering | CZ Biohub SF",
    author_email="michelle.khoo@czbiohub.org",
    license="MIT",
    packages=find_packages(),
<<<<<<< HEAD
    package_data={'stats_utils': ['stats_utils/data_files/*']},
=======
>>>>>>> parent of 5ec8a2d (Updated package installation settings)
    install_requires=[
        "numpy",
        "pandas",
        "black==23.1.0",
        "mypy==1.0.1",
        "mypy-extensions==1.0.0",
        "ruff==0.0.253",
    ],
    classifiers=["CZ Biohub :: Bioengineering"],
)
