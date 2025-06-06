#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: mwlkhoo
"""

from setuptools import setup, find_packages
from pathlib import Path


def readme():
    with open("README.md") as f:
        return f.read()


def get_data_files(parent_dir):
    parent_path = Path(parent_dir)
    child_paths = [path for path in parent_path.iterdir() if path.is_dir()]

    all_files = []
    for child_path in child_paths:
        files = [
            path.resolve().as_posix() for path in child_path.iterdir() if path.is_file()
        ]
        all_files.extend(files)

    return all_files


data_files = get_data_files("stats_utils/data_files")

setup(
    name="stats_utils",
    version="0.0.14",
    description="Statistics utilities for remoscope and corresponding paper",
    long_description=readme(),
    url="https://github.com/czbiohub-sf/remo-stats-utils",
    author="Bioengineering | CZ Biohub SF",
    author_email="michelle.khoo@czbiohub.org",
    license="BSD-3-Clause",
    packages=find_packages(),
    package_data={"stats_utils": data_files},
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
