import sys
from os.path import join, dirname, abspath

from setuptools import setup, find_packages

# this is a trick to get the version before the package is installed
directory = dirname(abspath(__file__))
sys.path.insert(0, join(directory, "bigg_models"))
version = "0.0.1"

setup(
    name="biggr_maps",
    version=version,
    author="Pascal Pieters",
    author_email="paspie@biosustain.dtu.dk",
    packages=find_packages(),
    package_data={
        "biggr_maps": []
    },
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.10",
    ],
    install_requires=[
        "escher",
    ],
)