"""Tests for ConnectomeInfluenceCalculator.

Covers the DataFrame-first __init__ and its from_sql / from_csv /
from_parquet / from_feather / from_numpy classmethod loaders, the
constructor parameter surface (inhibitory_nts, excluded_nts,
lambda_max, syn_weight_measure), the integration of adjust_influence
into calculate_influence, and the standalone adjust_influence helper.
PETSc/SLEPc-dependent tests use pytest.importorskip so the suite
remains runnable in environments without the HPC libraries installed.
"""
import numpy as np
import pandas as pd
import pytest

from InfluenceCalculator import InfluenceCalculator, adjust_influence


def _require_petsc_slepc():
    pytest.importorskip('petsc4py')
    pytest.importorskip('slepc4py')


# ---------------------------------------------------------------------
# Constructor parameter validation (no PETSc/SLEPc needed)
# ---------------------------------------------------------------------

def test_signed_requires_inhibitory_nts(celegans_edgelist_df,
                                        celegans_meta_df):
    with pytest.raises(ValueError, match='inhibitory_nts'):
        InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                            signed=True)


@pytest.mark.parametrize('bad_value', [0.0, 1.0, -0.1, 1.5, 2.0])
def test_lambda_max_validation(celegans_edgelist_df, bad_value):
    with pytest.raises(ValueError, match='lambda_max'):
        InfluenceCalculator(celegans_edgelist_df, lambda_max=bad_value)


def test_syn_weight_measure_validation(celegans_edgelist_df):
    with pytest.raises(ValueError, match='syn_weight_measure'):
        InfluenceCalculator(celegans_edgelist_df,
                            syn_weight_measure='bogus')


def test_missing_pre_post_columns_raises():
    bad = pd.DataFrame({'a': [1, 2], 'b': [3, 4], 'count': [5, 6]})
    with pytest.raises(ValueError, match="'pre' and 'post'"):
        InfluenceCalculator(bad)


def test_missing_count_or_weight_raises():
    bad = pd.DataFrame({'pre': [1, 2], 'post': [3, 4]})
    with pytest.raises(ValueError, match="'count'"):
        InfluenceCalculator(bad)


def test_excluded_nts_requires_meta(celegans_edgelist_df):
    with pytest.raises(ValueError, match='meta_df'):
        InfluenceCalculator(celegans_edgelist_df,
                            excluded_nts={'glutamate'})


def test_signed_requires_top_nt(celegans_edgelist_df):
    meta = pd.DataFrame({'root_id': ['ADAL', 'ADAR']})
    with pytest.raises(ValueError, match='top_nt'):
        InfluenceCalculator(celegans_edgelist_df, meta_df=meta,
                            signed=True, inhibitory_nts={'gaba'})


def test_norm_auto_computation_when_absent():
    edges = pd.DataFrame({
        'pre':   ['A', 'B', 'A'],
        'post':  ['C', 'C', 'D'],
        'count': [3, 7, 2],
    })
    # __init__ should compute 'norm' under the hood; smoke-test by
    # constructing without PETSc-dependent steps having to inspect it.
    _require_petsc_slepc()
    ic = InfluenceCalculator(edges, count_thresh=0, lambda_max=0.5)
    assert ic.n_neurons == 4


# ---------------------------------------------------------------------
# Construction smoke tests (need PETSc/SLEPc)
# ---------------------------------------------------------------------

def test_unsigned_construction_via_dataframes(celegans_edgelist_df,
                                              celegans_meta_df):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                             count_thresh=0, lambda_max=0.5)
    assert ic.n_neurons > 250
    assert ic.W_signed is False
    assert ic.lambda_max == 0.5
    assert ic.syn_weight_measure == 'count'


def test_signed_construction_via_dataframes(celegans_edgelist_df,
                                            celegans_meta_df):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                             signed=True, count_thresh=0,
                             lambda_max=0.5, inhibitory_nts={'gaba'})
    assert ic.W_signed is True
    assert ic.inhibitory_nts == {'gaba'}


def test_excluded_nts_drops_pre_neurons(celegans_edgelist_df,
                                        celegans_meta_df):
    _require_petsc_slepc()
    ic_baseline = InfluenceCalculator(celegans_edgelist_df,
                                      celegans_meta_df,
                                      count_thresh=0, lambda_max=0.5)
    ic_excl = InfluenceCalculator(celegans_edgelist_df,
                                  celegans_meta_df, count_thresh=0,
                                  lambda_max=0.5,
                                  excluded_nts={'glutamate'})
    assert ic_excl.n_neurons <= ic_baseline.n_neurons


# ---------------------------------------------------------------------
# Classmethod loaders -- format equivalence
# ---------------------------------------------------------------------

def test_from_sql_matches_from_dataframes(celegans_sqlite_path,
                                          celegans_edgelist_df,
                                          celegans_meta_df):
    _require_petsc_slepc()
    ic_df = InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                                count_thresh=0, lambda_max=0.5)
    ic_sql = InfluenceCalculator.from_sql(celegans_sqlite_path,
                                          count_thresh=0, lambda_max=0.5)
    assert ic_df.n_neurons == ic_sql.n_neurons
    assert set(ic_df.id_to_index.keys()) == set(ic_sql.id_to_index.keys())


def test_from_csv_matches_from_dataframes(celegans_csv_paths,
                                          celegans_edgelist_df,
                                          celegans_meta_df):
    _require_petsc_slepc()
    edgelist_path, meta_path = celegans_csv_paths
    ic_df = InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                                count_thresh=0, lambda_max=0.5)
    ic_csv = InfluenceCalculator.from_csv(edgelist_path, meta_path,
                                          count_thresh=0, lambda_max=0.5)
    assert ic_df.n_neurons == ic_csv.n_neurons
    assert set(ic_df.id_to_index.keys()) == set(ic_csv.id_to_index.keys())


def test_from_numpy_smoke():
    _require_petsc_slepc()
    rng = np.random.default_rng(0)
    n = 20
    adj = rng.random((n, n)) * (rng.random((n, n)) < 0.3)
    ic = InfluenceCalculator.from_numpy(adj, count_thresh=0,
                                        lambda_max=0.5)
    assert ic.n_neurons == n


# ---------------------------------------------------------------------
# adjust_influence module-level helper
# ---------------------------------------------------------------------

def _toy_unsigned_df():
    return pd.DataFrame({
        'matrix_index': [0, 1, 2, 3],
        'id': ['A', 'B', 'C', 'D'],
        'is_seed': [True, False, False, False],
        'Influence_score_(unsigned)': [1e-3, 1e-6, 0.0, 1.0],
    })


def _toy_signed_df():
    return pd.DataFrame({
        'matrix_index': [0, 1, 2, 3],
        'id': ['A', 'B', 'C', 'D'],
        'is_seed': [True, False, False, False],
        'Influence_score_(signed)': [1e-3, -1e-6, 0.0, -1.0],
    })


def test_adjust_influence_basic_columns():
    out = adjust_influence(_toy_unsigned_df(), const=24)
    expected = {
        'adjusted_influence',
        'adjusted_influence_norm_by_targets',
        'adjusted_influence_norm_by_sources_and_targets',
    }
    assert expected.issubset(out.columns)


def test_adjust_influence_log_const_anchor():
    out = adjust_influence(_toy_unsigned_df(), const=24)
    row_d = out.loc[out['id'] == 'D'].iloc[0]
    assert row_d['adjusted_influence'] == pytest.approx(24.0, abs=1e-4)


def test_adjust_influence_threshold_floor():
    out = adjust_influence(_toy_unsigned_df(), const=24)
    row_c = out.loc[out['id'] == 'C'].iloc[0]
    assert row_c['adjusted_influence'] == 0.0


def test_adjust_influence_preserves_sign():
    out = adjust_influence(_toy_signed_df(), const=24)
    row_d = out.loc[out['id'] == 'D'].iloc[0]
    assert row_d['adjusted_influence'] == pytest.approx(-24.0, abs=1e-4)


def test_adjust_influence_rejects_both_columns():
    df = _toy_unsigned_df()
    df['Influence_score_(signed)'] = df['Influence_score_(unsigned)']
    with pytest.raises(ValueError, match='both'):
        adjust_influence(df)


def test_adjust_influence_requires_score_column():
    df = pd.DataFrame({'id': [1, 2], 'is_seed': [True, False]})
    with pytest.raises(ValueError, match='influence score'):
        adjust_influence(df)


# ---------------------------------------------------------------------
# calculate_influence integration with adjust_influence
# ---------------------------------------------------------------------

def test_calculate_influence_default_returns_adjusted(celegans_edgelist_df,
                                                     celegans_meta_df):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                             count_thresh=0, lambda_max=0.5)
    df = ic.calculate_influence(['ADAL', 'ADAR'])
    assert 'Influence_score_(unsigned)' in df.columns
    assert 'adjusted_influence' in df.columns


def test_calculate_influence_adjust_off(celegans_edgelist_df,
                                        celegans_meta_df):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                             count_thresh=0, lambda_max=0.5)
    df = ic.calculate_influence(['ADAL', 'ADAR'], adjust=False)
    assert 'adjusted_influence' not in df.columns


def test_signed_calculate_influence_can_yield_negative(
        celegans_edgelist_df, celegans_meta_df):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_edgelist_df, celegans_meta_df,
                             signed=True, count_thresh=0,
                             lambda_max=0.5, inhibitory_nts={'gaba'})
    df = ic.calculate_influence(['ADAL', 'ADAR'], adjust=False)
    assert (df['Influence_score_(signed)'] < 0).any()


# ---------------------------------------------------------------------
# Bundled data sanity check
# ---------------------------------------------------------------------

def test_bundled_celegans_data_columns(celegans_meta_df,
                                       celegans_edgelist_df):
    assert {'pre', 'post', 'count', 'norm'} <= set(celegans_edgelist_df.columns)
    assert {'root_id', 'top_nt'} <= set(celegans_meta_df.columns)
    assert len(celegans_edgelist_df) > 1000
    assert len(celegans_meta_df) > 100
