
### Description
The Electronic Chirality Measure (ECM) is a quantity that allows to measure
the chirality of any (chiral) molecular system. It was proposed by Luca Bellarosa and Francesco
Zerbetto in 2003 [[1]](#1). In 2023 [[2]](#2) J. J. Aucar, A. Stroppa and G. A. Aucar have shown a novel, strong, 
and positive correlation between the energy difference of the total electronic energies of two enantiomers and ECM,
 supporting a subtle interplay between the weak forces acting within the nuclei of a given molecule and its chirality.


The Electronic Chirality Measure (ECM) is an end-to-end package implemented in Python 3.9 to 
measure the mentioned quantity. It also has some development interface
with the [PySCF](https://pyscf.org/) and [DIRAC](https://www.diracprogram.org) packages.

The ecm package uses simple input files in xyz format. The so called nearest
asymmetric structure is also needed. It can be easily obtained from (citar).
The package also allows the user to perform some simple plots to have a 
better understanding on the geometry of the systems under study.

The ECM is calculated computing the electronic wave function using a standard and powerful
quantum chemistry package such as [PySCF](https://pyscf.org/), the 
Python-based Simulations of Chemistry Framework. PySCF is an efficient 
platform for quantum chemistry calculations that can be used to simulate the
properties of molecules and crystals using mean-field and post-mean-field 
methods [[3]](#3).

The ecm main purpose is to be a user friendly package, easy to install, 
import, and run, available in most platforms and open-source. 
As a Python module, ecm can be integrated into any workflow. This package 
has code reviews with unit testing and continuous integration, code coverage
tools, and automatically keeps documentation up–to–date. 

The basic [example](docs/source/quickstart.rst) demonstrates some of these 
features.

This package allows to:

   1. Perform basic plots in order to analyze the systems under study.
   2. Import and export molecules from/to xyz and DIRAC format.
   3. Calculate CCM and ECM in a single structure.
   4. Calculate CCM and ECM in several structures at once and in the virtual mirroring path of each one.

Read more on [ECM publications](https://pubs.acs.org/doi/pdf/10.1021/ja028646%2B).

### Technical Documentation
Technical documents behind this project can be accessed [here](https://juanjoaucar.gitlab.io/pyECM).


### Requirements

First you should install the required python packages 

      - jupyterlab==3.1.13
      - jupyter==1.0.0
      - ipython[all]
      - attrs==22.2.0
      - numpy==1.21.2
      - pyscf==2.1.1
      - matplotlib==3.4.2
      - mendeleev==0.12.1

check the file `requirements.txt`. For developer, you should install `requirements_dev.txt`.

### Installation
pyECM is a **Python 3.9** package

1. Install virtual environment:

    ```python -m venv venv```

2. Activate virtual environment:

    ```source venv/bin/activate```

3. Install the packages:

      ```pip install pyECM```

4. Run pyECM (check some examples below)
    
5. For developer only, install dependencies:

      ```pip install -r requirements.txt -r requirements_dev.txt```

6. Run all test:

      ``tox==3.24.3``

### Usage

A detail workflow is provide into `workflow` directory. It has a list of Jupyter notebook with detail examples about pyECM tools and capabilities.

1. Getting starting with some basic plots.
      - `01.Plots.ipynb` 
2. Calculation of ECM in a simple molecule. 
      - `02.ECM_one-molecule.ipynb` 
3. Calculation of ECM in several molecules. 
      - `03.ECM_several-molecules.ipynb` 

## Roadmap
Some of the ideas to keep growing are:

* Obtain the nearest achiral structure within the code.

## Contributing
The easiest way to get help with the project is through the github project.

- GitHub:  git@github.com:juanjoaucar/pyECM.git

### Licence
GNU General Public License v3 (GLPv3)

### Authors and Acknowledgment
Main authors: Juan Jose Aucar (_juanaucar@gmail.com_)

Advisors: Gustavo A. Aucar and Alessandro Stroppa

### Project Status

Under development


Citing PyECM
------------
The following should be cited in publications utilizing the PyECM program package:

[A Relationship between the Molecular Parity-Violation Energy and the Electronic Chirality Measure](https://pubs.acs.org/doi/10.1021/acs.jpclett.3c03038),
J. J. Aucar, A. Stroppa, G. A. Aucar (2023),
*J. Phys. Chem. Lett.* doi:[10.1021/acs.jpclett.3c03038](https://pubs.acs.org/doi/10.1021/acs.jpclett.3c03038)

[PyECM23](https://doi.org/10.5281/zenodo.10149807)
Aucar, J. J. (2023),
*Zenodo*. doi: [10.5281/zenodo.10149807](https://doi.org/10.5281/zenodo.10149807)

[PySCF: the Python‐based simulations of chemistry framework](https://onlinelibrary.wiley.com/doi/abs/10.1002/wcms.1340),
Q. Sun, T. C. Berkelbach, N. S. Blunt, G. H. Booth, S. Guo, Z. Li, J. Liu,
J. McClain, E. R. Sayfutyarova, S. Sharma, S. Wouters, G. K.-L. Chan (2018),
*WIREs Comput. Mol. Sci.*, **8**: e1340. doi:[10.1002/wcms.1340](https://onlinelibrary.wiley.com/doi/abs/10.1002/wcms.1340)



Bug reports and feature requests
--------------------------------
Please submit tickets on the [issues](https://github.com/juanjoaucar/pyECM/issues) page.



---

### References
<div style=font-size:12px>
      <a id="1">[1]</a> 
      Bellarosa, L., & Zerbetto, F. (2003). Enantiomeric excesses and electronic chirality measure. Journal of the American Chemical Society, 125(7), 1975-1979.
<br>
</div>
<div style=font-size:12px>
      <a id="2">[2]</a> 
      J. J. Aucar, A. Stroppa, G. A. Aucar (2033). A Relationship between the Molecular Parity-Violation Energy and the Electronic Chirality Measure. Journal of Physical Chemistry Letters.
<br>
</div>
<div style=font-size:12px>
      <a id="3">[3]</a> 
      SUN, Qiming, et al. PySCF: the Python‐based simulations of chemistry framework. Wiley Interdisciplinary Reviews: Computational Molecular Science, 2018, vol. 8, no 1, p. e1340.
<br>
</div>

