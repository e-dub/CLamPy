<p align=center><img width="75%" src="figures/CLamPy.png"></p>

[![PyPi Version](https://img.shields.io/pypi/v/CLamPy.svg?style=flat-square)](https://pypi.org/project/CLamPy)
[![PyPI pyversions](https://img.shields.io/pypi/pyversions/CLamPy.svg?style=flat-square)](https://pypi.org/project/CLamPy/)
[![GitHub stars](https://img.shields.io/github/stars/e-dub/CLamPy.svg?style=flat-square&logo=github&label=Stars&logoColor=white)](https://github.com/e-dub/CLamPy)
[![PyPi downloads](https://img.shields.io/pypi/dm/CLamPy.svg?style=flat-square)](https://pypistats.org/packages/CLamPy)
[![Code style: blue](https://img.shields.io/badge/code%20style-blue-blue.svg)](https://blue.readthedocs.io/)

# CLamPy 

**Classical LAMinate theory for the lightweight design of structures and systems in PYthon**

**Klassische Laminattheorie für Leichtbau-Strukturen und -Systeme in Python**

**Teoria classica dei laminati per la costruzione leggera di strutture e sistemi in Python**


## Installation
### Prerequisites
Python 3 in addition to further packages. These necessary libraries can be installed via PIP:
```
pip install scipy
pip install numpy
pip install matplotlib
pip install tikzplotlib
pip install pandas
pip install openpyxl
```

### Install
```
python setup.py install
```

### PIP
You can also install CLamPy via PIP
```
pip install CLamPy
```

## Getting started
See Python scripts and Jupyter notebook in test.

## Assumptions
* Kirchhoff kinematic assumptions
  * Normals to the neutral plane remain normal after deformation
  * Normals remain straight after deformation
  * Thickness remains the same after deformation
* Ideal bonding
  *  No height of bonding between plies
  * Plies cannot slip relative to each other, i.e. no shear deformation
  * The strength between plies is infinite
* Geometry
  * Constant uniform thickness
  * Thin laminate, i.e. width and length > 10 times thickness
  * Small displacements and rotations, i.e. out-of-place displacement << thickness

## Furture work
- [ ] Add documentation including list of all nomenclature with meaning
- [ ] Different material parameters for each ply as array
- [ ] Material library
- [ ] Plot Young's modulus of each ply
- [ ] Top and Bottom characterisctics in an array for each ply
- [ ] Set warning if material properites have strange values
- [ ] Add thermal properties
- [ ] Add moisture properties

## Main source
* Jones, R. M. (2014) Mechanics of composite materials. Second Edition, Taylor and Francis, New York.

## Further references
* http://kmjespersen.sci-life.com/laminate-theory-example-using-python-notebook/
* https://github.com/joaopbernhardt/lamipy
* https://github.com/Eacaen/CLT-material-properties
* https://github.com/kemeen/material_mechanics
* http://kmjespersen.sci-life.com/laminate-theory-example-using-python-notebook/
* https://wstein.org/edu/2010/480b/projects/05-lamination_theory/
    
## Authors
E. M. Gioia & E. J. Wehrle

