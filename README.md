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

Another alternative would be to install both libraries using Conda. In this case, we highly recommend creating a virtual environment with a specific Python version:
```sh
conda create -n ic-venv-1 python=3.13.1
```
The packages can then be installed by executing the following:
```sh
conda install -c conda-forge petsc petsc4py
conda install -c conda-forge slepc slepc4py
```

## Description

This code computes the influence scores of a neuron or a group of neurons (seed) on all downstream neurons in the connectome based on a linear dynamical model of neural signal propagation: 

$$
\begin{aligned}
\tau \frac{d \boldsymbol{r}(t)}{dt} &= - \boldsymbol{r}(t) + \boldsymbol{W} \boldsymbol{r}(t) + \boldsymbol{s}(t) \\
&= \left( \boldsymbol{W} - \boldsymbol{I} \right)\boldsymbol{r}(t) + \boldsymbol{s}(t)
\end{aligned}
$$

where $\boldsymbol{r}$ is the vector of neural activity, $\boldsymbol{W}$ is the connectivity matrix, and $\boldsymbol{s}$ is the stimulus signal (applied to the seed neurons and remains constant throughout the simulation). The connectivity matrix is constructed by mapping the neuron IDs to matrix indices and arranging them such that the columns correspond to presynaptic neurons and rows correspond to postsynaptic neurons. The matrix is then filled with the number of synaptic connections that a presynaptic neuron projects onto a postsynaptic neuron.

To ensure stability of the dynamics, we rescale $\boldsymbol{W}$ such that:

$$
\begin{equation}
    \tilde{\boldsymbol{W}} = \frac{0.99}{\lambda_{max}} \boldsymbol{W}
\end{equation}
$$

where $\lambda_{max}$ is the largest eigenvalue of $\boldsymbol{W}$. The steady-state solution can thus be written as:

$$
\begin{equation}
    \boldsymbol{r}_{\infty} = -\boldsymbol{A}^{-1} \boldsymbol{s}
\end{equation}
$$

where $\boldsymbol{A} = \tilde{\boldsymbol{W}} - \boldsymbol{I}$. All matrix computations are performed using parallel computing libraries PETSc and SLEPc which adapt well to problems involving large, sparse matrices such as in our case and allow fast computation of the steady-state solution.

The influence of any seed is defined as the magnitude of neural activity at steady state. We wrote a Python code that outputs a Pandas dataframe with neuron IDs and the associated degree of influence a chosen seed exerts on them. The code also allows users to silence specific neurons throughout the simulation by setting appropriate entries of the connectivity matrix to zero, effectively cutting all synaptic connections from these neurons. This helps analyze the impact of any neuron along any pathway between seed and target neurons.

## Usage

To run a test example, start by importing the InfluenceCalculator package:

```python
from InfluenceCalculator import InfluenceCalculator
```
Then instantiate a class object 'ic' using the filepath to the BANC dataset (should be 'sqlite' file):
```python
# Build InfluenceCalculator object
ic = InfluenceCalculator('BANC_dataset.sqlite')
```

Let us now, define the seed group as all 'olfactory' neurons and calculate the influence of this seed on all downstream neurons while making sure to inhibit all non-seed sensory neurons:

```python
# Define seed category (depending on how neurons are labelled in metadata)
meta_column = 'seed_01'
seed_category = 'olfactory'

# Get seed neuron ids
seed_ids = ic.meta[ic.meta[meta_column] == seed_category].root_id 

# Get neuron ids to inhibit (sensory neurons in this case)
silenced_neurons = ic.meta[
    ic.meta['super_class'].isin(['sensory',
                                    'ascending_sensory'])].root_id

# Calculate influence scores and store them in a Pandas dataframe
influence_df = ic.calculate_influence(seed_ids, silenced_neurons)
```

Executing the previous script outputs a dataframe with a column of neurons IDs, a column of Boolean entries indicating whether the corresponding neuron is part of the seed group or not and, a column of influence scores relative to the seed neurons. Note that even though we selected all sensory neurons in the inhibition group, the 'calculate_influence' routine excludes the seed neurons from being inhibited if they intitially belonged to that group.


## Contributing

TODO

## Licence

TODO
