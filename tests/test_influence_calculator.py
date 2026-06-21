"""Tests for ConnectomeInfluenceCalculator.

Covers the constructor parameter surface introduced in the inhibitory_nts
/ excluded_nts / lambda_max / syn_weight_measure PR, the integration of
adjust_influence into calculate_influence, and the InfluenceCalculator
.adjust_influence @staticmethod itself.  PETSc/SLEPc-dependent tests
use pytest.importorskip so the suite remains runnable in environments
without the HPC libraries installed.
"""
import numpy as np
import pandas as pd
import pytest

from InfluenceCalculator import InfluenceCalculator


# ---------------------------------------------------------------------
# Constructor parameter validation (no PETSc/SLEPc needed)
# ---------------------------------------------------------------------

def test_signed_requires_inhibitory_nts(celegans_sqlite_path):
    with pytest.raises(ValueError, match='inhibitory_nts'):
        InfluenceCalculator(celegans_sqlite_path, signed=True)


@pytest.mark.parametrize('bad_value', [0.0, 1.0, -0.1, 1.5, 2.0])
def test_lambda_max_validation(celegans_sqlite_path, bad_value):
    with pytest.raises(ValueError, match='lambda_max'):
        InfluenceCalculator(celegans_sqlite_path, lambda_max=bad_value)


def test_syn_weight_measure_validation(celegans_sqlite_path):
    with pytest.raises(ValueError, match='syn_weight_measure'):
        InfluenceCalculator(celegans_sqlite_path,
                            syn_weight_measure='bogus')


# ---------------------------------------------------------------------
# Construction smoke tests (need PETSc/SLEPc)
# ---------------------------------------------------------------------

def _require_petsc_slepc():
    pytest.importorskip('petsc4py')
    pytest.importorskip('slepc4py')


def test_unsigned_construction_smoke(celegans_sqlite_path):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_sqlite_path, count_thresh=0,
                             lambda_max=0.5)
    assert ic.n_neurons > 250
    assert ic.W_signed is False
    assert ic.lambda_max == 0.5
    assert ic.syn_weight_measure == 'count'


def test_signed_construction_smoke(celegans_sqlite_path):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_sqlite_path, signed=True,
                             count_thresh=0, lambda_max=0.5,
                             inhibitory_nts={'gaba'})
    assert ic.W_signed is True
    assert ic.inhibitory_nts == {'gaba'}


def test_norm_syn_weight_measure_construction_smoke(celegans_sqlite_path):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_sqlite_path, count_thresh=0,
                             lambda_max=0.5, syn_weight_measure='norm')
    assert ic.syn_weight_measure == 'norm'


def test_excluded_nts_drops_pre_neurons(celegans_sqlite_path):
    _require_petsc_slepc()
    ic_baseline = InfluenceCalculator(celegans_sqlite_path,
                                      count_thresh=0, lambda_max=0.5)
    ic_excl = InfluenceCalculator(celegans_sqlite_path, count_thresh=0,
                                  lambda_max=0.5,
                                  excluded_nts={'glutamate'})
    # Excluding glutamate drops every edge whose pre-neuron is a
    # glutamate cell, so any neuron that only ever appeared as a
    # glutamatergic pre disappears from the W mapping entirely.
    assert ic_excl.n_neurons <= ic_baseline.n_neurons


# ---------------------------------------------------------------------
# InfluenceCalculator.adjust_influence @staticmethod
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
    out = InfluenceCalculator.adjust_influence(_toy_unsigned_df(), const=24)
    expected = {
        'adjusted_influence',
        'adjusted_influence_norm_by_targets',
        'adjusted_influence_norm_by_sources_and_targets',
    }
    assert expected.issubset(out.columns)


def test_adjust_influence_log_const_anchor():
    out = InfluenceCalculator.adjust_influence(_toy_unsigned_df(), const=24)
    # Strongest influence (1.0) gives log(1.0) + 24 = 24
    row_d = out.loc[out['id'] == 'D'].iloc[0]
    assert row_d['adjusted_influence'] == pytest.approx(24.0, abs=1e-4)


def test_adjust_influence_threshold_floor():
    out = InfluenceCalculator.adjust_influence(_toy_unsigned_df(), const=24)
    # 0.0 raw influence floors to exp(-24); log(exp(-24)) + 24 = 0
    row_c = out.loc[out['id'] == 'C'].iloc[0]
    assert row_c['adjusted_influence'] == 0.0


def test_adjust_influence_preserves_sign():
    out = InfluenceCalculator.adjust_influence(_toy_signed_df(), const=24)
    row_d = out.loc[out['id'] == 'D'].iloc[0]
    assert row_d['adjusted_influence'] == pytest.approx(-24.0, abs=1e-4)


def test_adjust_influence_rejects_both_columns():
    df = _toy_unsigned_df()
    df['Influence_score_(signed)'] = df['Influence_score_(unsigned)']
    with pytest.raises(ValueError, match='both'):
        InfluenceCalculator.adjust_influence(df)


def test_adjust_influence_requires_score_column():
    df = pd.DataFrame({'id': [1, 2], 'is_seed': [True, False]})
    with pytest.raises(ValueError, match='influence score'):
        InfluenceCalculator.adjust_influence(df)


# ---------------------------------------------------------------------
# calculate_influence integration with adjust_influence
# ---------------------------------------------------------------------

def test_calculate_influence_default_returns_adjusted(celegans_sqlite_path):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_sqlite_path, count_thresh=0,
                             lambda_max=0.5)
    df = ic.calculate_influence(['ADAL', 'ADAR'])
    assert 'Influence_score_(unsigned)' in df.columns
    assert 'adjusted_influence' in df.columns
    assert 'adjusted_influence_norm_by_targets' in df.columns
    assert 'adjusted_influence_norm_by_sources_and_targets' in df.columns


def test_calculate_influence_adjust_off(celegans_sqlite_path):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_sqlite_path, count_thresh=0,
                             lambda_max=0.5)
    df = ic.calculate_influence(['ADAL', 'ADAR'], adjust=False)
    assert 'Influence_score_(unsigned)' in df.columns
    assert 'adjusted_influence' not in df.columns


def test_signed_calculate_influence_can_yield_negative(celegans_sqlite_path):
    _require_petsc_slepc()
    ic = InfluenceCalculator(celegans_sqlite_path, signed=True,
                             count_thresh=0, lambda_max=0.5,
                             inhibitory_nts={'gaba'})
    df = ic.calculate_influence(['ADAL', 'ADAR'], adjust=False)
    raw = df['Influence_score_(signed)']
    # With gaba-mediated inhibition active some downstream targets
    # should accumulate net-negative influence.
    assert (raw < 0).any()


# ---------------------------------------------------------------------
# Bundled data sanity check
# ---------------------------------------------------------------------

def test_bundled_celegans_data_columns(celegans_meta_df, celegans_edgelist_df):
    assert {'pre', 'post', 'count', 'norm'} <= set(celegans_edgelist_df.columns)
    assert {'root_id', 'top_nt'} <= set(celegans_meta_df.columns)
    assert len(celegans_edgelist_df) > 1000
    assert len(celegans_meta_df) > 100
