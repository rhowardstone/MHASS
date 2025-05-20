from setuptools import setup, find_packages
import os

# Read the contents of README.md
with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="mhass",
    version="0.1.0",
    author="Rye Howard-Stone",
    author_email="Rye.Howard-Stone@UConn.edu",
    description="A complete pipeline for simulating PacBio ccs amplicon data",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/rhowardstone/MHASS",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "tqdm>=4.50.0",
    ],
    entry_points={
        "console_scripts": [
            "mhass=mhass.main:main",
        ],
    },
    package_data={
        "mhass": [
            "resources/*",
            "resources/pbsim3/src/pbsim",
            "resources/pbsim3/data/QSHMM-RSII.model",
            "get_counts.R",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    python_requires=">=3.6",
)
