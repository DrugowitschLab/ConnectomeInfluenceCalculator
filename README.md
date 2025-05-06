\usepackage{bm} 

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

This code computes the influence scores of a neuron or a group of neurons (seed) on all downstream neurons in the connectome based on a linear dynamical model of neural signal propagation: 

$$
\begin{aligned}
\tau \frac{d \boldsymbol{r}(t)}{dt} &= - \boldsymbol{r}(t) + \boldsymbol{W} \boldsymbol{r}(t) + \boldsymbol{s}(t) \\
&= \left( \boldsymbol{W} - \boldsymbol{I} \right)\boldsymbol{r}(t) + \boldsymbol{s}(t)
\end{aligned}
$$

where $\boldsymbol{r}$ is the vector of neural activity, $\boldsymbol{W}$ is the connectivity matrix, and $\boldsymbol{s}$ is the stimulus signal (applied to the seed neurons and remains constant throughout the simulation). The connectivity matrix is constructed by mapping the neuron IDs to matrix indices and arranging them such that the columns correspond to presynaptic neurons and rows correspond to postsynaptic neurons. The matrix is then filled with the number of synaptic connections that a presynaptic neuron projects onto a postsynaptic neuron.


## Contributing

TODO

## Licence

TODO
