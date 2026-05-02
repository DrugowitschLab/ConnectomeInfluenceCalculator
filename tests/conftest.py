"""Test fixtures for ConnectomeInfluenceCalculator.

The DataFrame-first __init__ added in this PR lets the bundled
C. elegans CSVs feed every test path directly.  We expose them four
ways: as DataFrames, as filesystem paths to the bundled CSVs, and as a
session-scoped temporary SQLite database for exercising the from_sql
classmethod.
"""
import sqlite3
from importlib.resources import as_file, files

import pytest

from InfluenceCalculator.data import celegans_edgelist, celegans_meta


@pytest.fixture
def celegans_meta_df():
    return celegans_meta()


@pytest.fixture
def celegans_edgelist_df():
    return celegans_edgelist()


@pytest.fixture(scope='session')
def celegans_csv_paths(tmp_path_factory):
    """Materialise the bundled CSVs as filesystem paths so from_csv can
    consume them.  importlib.resources.files() may return a zip-internal
    path for installed wheels; as_file() copies to disk if needed.
    """
    pkg = files('InfluenceCalculator.data')
    paths = tmp_path_factory.mktemp('celegans_csv')
    edgelist_path = paths / 'celegans_edgelist.csv'
    meta_path = paths / 'celegans_meta.csv'
    with as_file(pkg / 'celegans_edgelist.csv') as src:
        edgelist_path.write_bytes(src.read_bytes())
    with as_file(pkg / 'celegans_meta.csv') as src:
        meta_path.write_bytes(src.read_bytes())
    return str(edgelist_path), str(meta_path)


@pytest.fixture(scope='session')
def celegans_sqlite_path(tmp_path_factory):
    """Build a temporary SQLite database from the bundled C. elegans
    edge list and metadata, matching the schema from_sql expects.
    """
    edges = celegans_edgelist()
    meta = celegans_meta()

    db_path = tmp_path_factory.mktemp('celegans_db') / 'celegans.sqlite'
    conn = sqlite3.connect(db_path)
    meta.to_sql('meta', conn, index=False)
    edges.to_sql('edgelist_simple', conn, index=False)
    conn.close()
    return str(db_path)
