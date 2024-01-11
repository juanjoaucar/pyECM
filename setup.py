import os
import pathlib

from setuptools import setup

# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))


REQUIREMENTS = [
	"jupyterlab==3.1.13",
	"jupyter==1.0.0",
	"ipython[all]",
	"attrs==22.2.0",
	"numpy==1.22.0",
	"pyscf==2.1.1",
	"matplotlib==3.4.2",
	"mendeleev==0.12.1",
    "sphinxcontrib-bibtex==2.5.0",
    "scipy==1.10.1",
    "openpyxl==3.1.2"
]

with open(PATH / "pyECM" / "__init__.py") as fp:
    for line in fp.readlines():
        if line.startswith("__version__ = "):
            VERSION = line.split("=", 1)[-1].replace('"', "").strip()
            break

with open("README.md", "r") as readme:
    LONG_DESCRIPTION = readme.read()


# =============================================================================
# FUNCTIONS
# =============================================================================

setup(
    name="pyECM",
    version="0.1.1",
    author="""
    Juan José Aucar
    """,
    author_email="""
    juanaucar@gmail.com, 
    """,
    packages=["pyECM"],
    install_requires=REQUIREMENTS,
    license="The GPLv3 License",
    description="Electronic Chirality Measure",
    long_description_content_type="text/markdown",
    long_description=LONG_DESCRIPTION,
    url="https://github.com/juanjoaucar/pyECM",
    keywords=[
        "Electronic Chirality Measure",
        "Continous Chirality Measure",
        "Relativity",
    ],
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Education",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: Implementation :: CPython",
        "Topic :: Scientific/Engineering",
    ],
    include_package_data=True,
)
