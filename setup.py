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
    install_requires=[
        "numpy",
        "pandas",
    ],
    classifiers=["CZ Biohub :: Bioengineering"],
)

