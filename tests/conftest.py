"""Test fixtures for ConnectomeInfluenceCalculator.

The InfluenceCalculator constructor in this PR series still expects an
SQLite file (DataFrame / CSV constructors arrive in a later PR), so the
session-scoped fixture builds a temporary SQLite database from the
bundled C. elegans CSVs to exercise the SQLite path against a real
connectome.
"""
import sqlite3

import pytest

from InfluenceCalculator.data import celegans_edgelist, celegans_meta


@pytest.fixture(scope='session')
def celegans_sqlite_path(tmp_path_factory):
    """Build a temporary SQLite database from the bundled C. elegans
    edge list and metadata, matching the schema InfluenceCalculator
    expects (a 'meta' table and an 'edgelist_simple' table).
    """
    edges = celegans_edgelist()
    meta = celegans_meta()

    db_path = tmp_path_factory.mktemp('celegans_db') / 'celegans.sqlite'
    conn = sqlite3.connect(db_path)
    meta.to_sql('meta', conn, index=False)
    edges.to_sql('edgelist_simple', conn, index=False)
    conn.close()
    return str(db_path)


@pytest.fixture
def celegans_meta_df():
    return celegans_meta()


@pytest.fixture
def celegans_edgelist_df():
    return celegans_edgelist()
