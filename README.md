# ConnectomeInfluenceCalculator

This code computes the influence scores of a neuron or a group of neurons (seed) on all downstream neurons in the connectome based on a linear dynamical model of neural signal propagation.

## Install

Download the repository, using either the [current main branch](https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator/archive/refs/heads/main.zip) or one of the [past releases](https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator/releases). Then, either directly install the zipped version using
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

Alteratively, both libraries and their Python wrappers can be installed using `conda`. In this case, we highly recommend creating a virtual environment with a specific Python version (version `3.13.1` worked for us):
```sh
conda create -n ic-venv python=3.13.1
```
After activation of the virtual environment, the packages can be installed by executing the following commands:
```sh
conda install -c conda-forge petsc petsc4py
conda install -c conda-forge slepc slepc4py
```

## Description

This code computes the influence scores of a neuron or a group of neurons, as specified through the seed vector $\boldsymbol{s}$, on all downstream neurons in the connectome based on a linear dynamical model of neural signal propagation: 

$$
\begin{aligned}
\tau \frac{d \boldsymbol{r}(t)}{dt} &= - \boldsymbol{r}(t) + \boldsymbol{W} \boldsymbol{r}(t) + \boldsymbol{s}(t) \\
&= \left( \boldsymbol{W} - \boldsymbol{I} \right)\boldsymbol{r}(t) + \boldsymbol{s}(t)
\end{aligned}
$$

where $\boldsymbol{r}$ is the vector of neural activity, $\boldsymbol{W}$ is the connectivity matrix, and $\boldsymbol{s}$ is the simulated neural stimulation (applied to the seed neurons and remains constant throughout the simulation). The connectivity matrix is constructed by mapping the neuron IDs to matrix indices and arranging them such that the columns correspond to presynaptic neurons and rows correspond to postsynaptic neurons. The matrix is then filled with the number of synaptic connections that a presynaptic neuron projects onto a postsynaptic neuron.

To ensure stable neural dynamics, we rescale $\boldsymbol{W}$ such that:

$$
\begin{equation}
    \tilde{\boldsymbol{W}} = \frac{0.99}{\lambda_{max}} \boldsymbol{W}
\end{equation}
$$

where $\lambda_{max}$ is the largest real eigenvalue of $\boldsymbol{W}$. The steady-state solution can thus be written as:

$$
\begin{equation}
    \boldsymbol{r}_{\infty} = - \left( \tilde{\boldsymbol{W}} - \boldsymbol{I} \right)^{-1} \boldsymbol{s} .
\end{equation}
$$

All matrix computations are performed using parallel computing libraries PETSc and SLEPc which adapt well to problems involving large, sparse matrices, like neural connectivity matrices, and allow fast computation of the steady-state solution.

The influence of any seed is defined as the magnitude of neural activity at steady state, $\boldsymbol{r}_{\infty}$.
Our Python code supports creating a Pandas dataframe with neuron IDs and the associated degree of influence a chosen seed exerts on them. The code also allows users to silence specific neurons throughout the simulation by setting appropriate entries of the connectivity matrix to zero, effectively cutting all synaptic connections from these neurons. This helps analyze the impact of any neuron along any pathway between seed and target neurons.

## Usage

To run a test example, start by importing the InfluenceCalculator package:

```python
from InfluenceCalculator import InfluenceCalculator
```
Then instantiate a class object `ic` using the filepath to the BANC connectome dataset (should be an SQLite file):
```python
# Build InfluenceCalculator object
ic = InfluenceCalculator('BANC_dataset.sqlite')
```

By default, the programme simulates neural signal propagation based on an unsigned version of the connectivity matrix. However, it is possible to use a signed version whereby synaptic weights of inhibitory neurons are assigned negative values. Moreover, users can specify a minimum threshold count for the number of postsynaptic connections that will be considered in the analysis (default is `5`). To do so, run the following command, instead:
```python
# Build InfluenceCalculator object
ic = InfluenceCalculator('BANC_dataset.sqlite', signed=True, count_thresh=5)
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

Executing this script returns a dataframe with a column of neurons IDs, a column of Boolean entries indicating whether the corresponding neuron is part of the seed group or not, and a column of influence scores relative to the seed neurons.
Note that even though we selected all sensory neurons to be inhibited, the `calculate_influence` method ensures that no seed neurons are being inhibited.


## Contributing

We love to hear your input, and want to make contributing to this project as easy and transparent as possible, whether it's:

- Reporting a bug
- Discussing the current state of the code
- Submitting a fix
- Proposing new features
- Becoming a maintainer

When contributing to this repository, please first discuss the change you wish to make via issue, email, or any other method with the owners of this repository before making a change.

## We Develop with Github

We use Github to host code, to track [issues and feature requests](https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator/issues), as well as accept [pull requests](https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator/pulls). 

## Any contributions you make will be under the New BSD License

In short, when you submit code changes, your submissions are understood to be under the same [New BSD License](https://choosealicense.com/licenses/bsd-3-clause/) that covers the project. Feel free to contact the maintainers if that's a concern.

## Report bugs using Github's [issues](https://github.com/DrugowitschLab/ConnectomeInfluenceCalculator/issues)

We use GitHub issues to track public bugs. Report a bug by opening a new issue; it's that easy!

If you file a bug report, please make sure to include

- A quick summary and/or background
- Steps to reproduce
  - Be specific!
  - Give minimal working example code that reproduces the bug, if you can.
- What you expected would happen
- What actually happens
- Notes (possibly including why you think this might be happening, or stuff you tried that didn't work)

# Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project maintainers at jan_drugowitsch@hms.harvard.edu 
or zaki_ajabi@hms.harvard.edu. All complaints will be reviewed and investigated 
and will result in a response that is deemed necessary and appropriate to the 
circumstances. The project maintainer is obligated to maintain confidentiality 
with regard to the reporter of an incident. Further details of specific enforcement 
policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
