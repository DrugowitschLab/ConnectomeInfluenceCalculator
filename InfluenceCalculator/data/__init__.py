"""Bundled connectome datasets for ConnectomeInfluenceCalculator.

C. elegans
----------
celegans_edgelist() and celegans_meta() return the C. elegans
hermaphrodite chemical connectome (300 cells, 3,539 edges, 20,672
synapses) as pandas DataFrames.

The CSVs were taken from the OpenWorm project distribution
(https://openworm.org/, accessed February 2026), which aggregates the
original electron-microscopy reconstructions of White et al. 1986 and
Cook et al. 2019.  Neurotransmitter and cell-class annotations follow
WormAtlas (https://www.wormatlas.org/) and CenGen (Taylor et al. 2021,
https://cengen.shinyapps.io/CengenApp/).

If you publish results computed from these CSVs, please cite the
OpenWorm aggregation alongside the underlying connectome papers:

    @article{szigeti2014openworm,
      title   = {{OpenWorm}: an open-science approach to modeling
                 Caenorhabditis elegans},
      author  = {Szigeti, Bal{\\'a}zs and Gleeson, Padraig and
                 Vella, Michael and others},
      journal = {Frontiers in Computational Neuroscience},
      volume  = {8}, pages = {137}, year = 2014,
      doi     = {10.3389/fncom.2014.00137}}

    @article{white1986structure,
      title   = {The structure of the nervous system of the nematode
                 Caenorhabditis elegans},
      author  = {White, J. G. and Southgate, E. and
                 Thomson, J. N. and Brenner, S.},
      journal = {Philosophical Transactions of the Royal Society B},
      volume  = {314}, number = {1165}, pages = {1--340}, year = 1986,
      doi     = {10.1098/rstb.1986.0056}}

    @article{cook2019whole,
      title   = {Whole-animal connectomes of both Caenorhabditis
                 elegans sexes},
      author  = {Cook, Steven J. and Jarrell, Travis A. and
                 Brittin, Christopher A. and others},
      journal = {Nature}, volume = {571}, number = {7763},
      pages   = {63--71}, year = 2019,
      doi     = {10.1038/s41586-019-1352-7}}
"""
import pandas as pd
from importlib.resources import files


def celegans_edgelist():
    """Load the bundled C. elegans synaptic edge list as a pandas
    DataFrame with columns 'pre', 'post', 'count', and 'norm'.
    """
    return pd.read_csv(files(__package__) / "celegans_edgelist.csv")


def celegans_meta():
    """Load the bundled C. elegans neuron metadata as a pandas
    DataFrame with columns 'root_id', 'top_nt', 'super_class',
    'neuron_class', and 'body_part'.
    """
    return pd.read_csv(files(__package__) / "celegans_meta.csv")


__all__ = ["celegans_edgelist", "celegans_meta"]
