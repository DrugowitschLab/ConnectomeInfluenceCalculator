import pytest
from importlib.resources import files, as_file

from InfluenceCalculator.data import celegans_edgelist, celegans_meta


@pytest.fixture(scope="session")
def edgelist_path():
    with as_file(files("InfluenceCalculator.data")
                 / "celegans_edgelist.csv") as p:
        yield p


@pytest.fixture(scope="session")
def meta_path():
    with as_file(files("InfluenceCalculator.data")
                 / "celegans_meta.csv") as p:
        yield p


@pytest.fixture(scope="session")
def edgelist_df():
    return celegans_edgelist()


@pytest.fixture(scope="session")
def meta_df():
    return celegans_meta()
