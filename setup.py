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

def package_files(directory):
    paths = []
    for (path, directories, filenames) in os.walk(directory):
        for filename in filenames:
            paths.append(os.path.join('..', path, filename))
    return paths

extra_files = package_files('stats_utils/data_files/frightful-wendigo-1931')
print(extra_files)

setup(
    name="stats_utils",
    version="0.0.3",
    description="Statistics utilities for remoscope and corresponding paper",
    long_description=readme(),
    url="https://github.com/czbiohub-sf/remo-stats-utils",
    author="Bioengineering | CZ Biohub SF",
    author_email="michelle.khoo@czbiohub.org",
    license="MIT",
    packages=find_packages(),
    package_data={'stats_utils': extra_files},
    include_package_data=True,
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
