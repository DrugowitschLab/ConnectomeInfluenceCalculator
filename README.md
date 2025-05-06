# ConnectomeInfluenceCalculator

## Install

Download the repository, using either the [current main branch](https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator/archive/refs/heads/main.zip) or one of the past versions. Then, either directly install the zipped version using
```sh
python3 -m pip install ConnectomeInfluenceCalculator-main.zip
```
or unzip the files into a folder and run
```sh
python3 -m pip install .
```
in that folder.

### Troubleshooting PETSc and SLEPc installations

The package relies on the the [PETSc](https://petsc.org/) and [SLEPc](https://slepc.upv.es) libraries to perform sparse matrix computations. These libraries consist of the core libraries, and associated Python wrappers. The core libraries might not install correctly through `pip`. If this happens, a possible work-around is to install them through other means (e.g., using Homebrew on OS X), and then tell the Python wrappers `petsc4py` and `slepc4py` where to find them by calling
```sh
export PETSC_DIR = /path/to/PETSc/installation
export SLEPC_DIR = /path/to/SLEPc/installation
```
before running the above `pip` commands. Please make sure that installed core libraries have the same version numbers as the Python wrappers that will be installed.

## Usage

We developed a method to compute the influence scores of a neuron or a group of neurons (seed) on all downstream neurons in the connectome based on a linear dynamical model of neural signal propagation: 
$$
x = \frac{-b \pm \sqrt{b^2 - 4ac}}{2a}
$$

## Contributing

TODO

## Licence

TODO
