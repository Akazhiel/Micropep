"""
Micropeptide obtention pipeline.
@author: Jonatan Gonzalez Rodriguez <jonatan.gonzalez.r@outlook.com>
"""

import os
import io
import glob
from setuptools import setup, find_packages

# Get the long description from the relevant file
here = os.path.abspath(os.path.dirname(__file__))
with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
    long_description = f.read()

setup(
    name="micropep",
    version="0.1.0",
    description=__doc__.split("\n", 1)[0],
    long_description=long_description,
    keywords="micropeptide obtention pipeline",
    author="Jonatan Gonzalez Rodriguez",
    author_email="jonatan.gonzalez.r@outlook.com",
    license="MIT",
    packages=find_packages(),
    include_package_data=False,
    package_data={"": ["RELEASE-VERSION"]},
    zip_safe=False,
    install_requires=[
        "setuptools",
        "numpy",
        "pandas",
        "BioPython",
        "pyranges",
    ],
    # test_suite = 'tests',
    scripts=glob.glob("*.py"),
    classifiers=[
        "Intended Audience :: Science/Research",
        "Topic :: Software Development",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: MIT:: Copyright Jonatan Gonzalez Rodriguez",
        "Programming Language :: Python :: 3.6",
        "Environment :: Console",
    ],
)