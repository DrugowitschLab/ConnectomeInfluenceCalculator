import pytest
import numpy as np
import pandas as pd

from InfluenceCalculator import InfluenceCalculator, adjust_influence
from InfluenceCalculator.InfluenceCalculator import _validate_and_prepare_edgelist


# ---------------------------------------------------------------------------
# 1. Format equivalence (requires PETSc)
# ---------------------------------------------------------------------------

def test_format_equivalence(edgelist_path, meta_path, edgelist_df, meta_df):
    """from_dataframes and from_csv produce identical internal state."""
    pytest.importorskip("petsc4py", reason="PETSc not installed")

    ic_df = InfluenceCalculator.from_dataframes(edgelist_df.copy(), meta_df)
    ic_csv = InfluenceCalculator.from_csv(str(edgelist_path), str(meta_path))

    # Build adjacency matrix for from_numpy (integer neuron indices)
    elist, _ = _validate_and_prepare_edgelist(edgelist_df.copy(), None, False, 5)
    unique_ids = pd.unique(
        np.hstack([elist["post"].values, elist["pre"].values])
    )
    n = len(unique_ids)
    id_to_idx = {nid: i for i, nid in enumerate(unique_ids)}
    W_dense = np.zeros((n, n))
    pre_idx = np.array([id_to_idx[p] for p in elist["pre"].values])
    post_idx = np.array([id_to_idx[p] for p in elist["post"].values])
    np.add.at(W_dense, (post_idx, pre_idx), elist["norm"].values)
    ic_np = InfluenceCalculator.from_numpy(W_dense)

    # All three must agree on neuron count and W dimensions
    assert ic_df.n_neurons == ic_csv.n_neurons == ic_np.n_neurons
    assert ic_df.W.getSize() == ic_csv.W.getSize() == ic_np.W.getSize()
    # from_dataframes and from_csv must share the same neuron ID universe
    assert set(ic_df.id_to_index.keys()) == set(ic_csv.id_to_index.keys())


# ---------------------------------------------------------------------------
# 2. adjust_influence basics
# ---------------------------------------------------------------------------

def test_adjust_influence_basics():
    """Output columns exist and computed values match expected log + const."""
    df = pd.DataFrame(
        {
            "id": ["A", "B", "C"],
            "is_seed": [True, False, False],
            "Influence_score_(unsigned)": [1.0, 0.5, 0.1],
        }
    )
    const = 24
    result = adjust_influence(df, const=const)

    expected_cols = [
        "adjusted_influence",
        "adjusted_influence_norm_by_targets",
        "adjusted_influence_norm_by_sources_and_targets",
    ]
    for col in expected_cols:
        assert col in result.columns, f"Missing column: {col}"

    # Without a 'target' column each neuron is its own group, so adjusted
    # influence for neuron A = log(1.0) + 24 = 24.0
    a_val = result.loc[result["id"] == "A", "adjusted_influence"].iloc[0]
    assert abs(a_val - 24.0) < 1e-4

    # norm_by_targets divides by n_targets (3 neurons) before log
    # For A: log(1.0 / 3) + 24 = log(1/3) + 24
    expected_nbt = np.log(1.0 / 3) + const
    a_nbt = result.loc[result["id"] == "A",
                        "adjusted_influence_norm_by_targets"].iloc[0]
    assert abs(a_nbt - expected_nbt) < 1e-4


# ---------------------------------------------------------------------------
# 3. adjust_influence threshold floor
# ---------------------------------------------------------------------------

def test_adjust_influence_threshold_floor():
    """Scores near zero must produce finite, non-negative adjusted values."""
    df = pd.DataFrame(
        {
            "id": ["A", "B", "C"],
            "is_seed": [True, False, False],
            "Influence_score_(unsigned)": [1e-200, 0.0, np.finfo(float).tiny],
        }
    )
    result = adjust_influence(df, const=24)

    for col in [
        "adjusted_influence",
        "adjusted_influence_norm_by_targets",
        "adjusted_influence_norm_by_sources_and_targets",
    ]:
        assert np.isfinite(result[col]).all(), f"{col} contains non-finite values"
        assert (result[col] >= 0).all(), f"{col} contains negative values"


# ---------------------------------------------------------------------------
# 3b. adjust_influence sign preservation (signed input)
# ---------------------------------------------------------------------------

def test_adjust_influence_preserves_sign():
    """Negative raw influence yields negative adjusted_influence; the
    floor still clips magnitudes below exp(-const) toward 0 in either
    sign."""
    df = pd.DataFrame(
        {
            "id": ["A", "B", "C", "D"],
            "is_seed": [True, False, False, False],
            "Influence_score_(signed)": [1.0, -0.5, 1e-200, -1e-200],
        }
    )
    const = 12
    result = adjust_influence(df, const=const).set_index("id")

    # A: log(1) + 12 = 12
    assert abs(result.loc["A", "adjusted_influence"] - 12.0) < 1e-4
    # B: -(log(0.5) + 12) ≈ -11.307
    expected_b = -(np.log(0.5) + const)
    assert abs(result.loc["B", "adjusted_influence"] - expected_b) < 1e-4
    # C and D are below the floor and should clip to ~0 regardless of sign
    assert abs(result.loc["C", "adjusted_influence"]) < 1e-4
    assert abs(result.loc["D", "adjusted_influence"]) < 1e-4


# ---------------------------------------------------------------------------
# 4. Input validation — missing required columns
# ---------------------------------------------------------------------------

def test_input_validation_missing_columns():
    """Missing 'pre' or 'post' in edgelist raises ValueError."""
    df_no_pre = pd.DataFrame({"post": [1, 2], "count": [5, 6]})
    with pytest.raises(ValueError, match="'pre' and 'post'"):
        InfluenceCalculator.from_dataframes(df_no_pre)

    df_no_post = pd.DataFrame({"pre": [1, 2], "count": [5, 6]})
    with pytest.raises(ValueError, match="'pre' and 'post'"):
        InfluenceCalculator.from_dataframes(df_no_post)


# ---------------------------------------------------------------------------
# 5. Input validation — signed=True without top_nt in meta
# ---------------------------------------------------------------------------

def test_input_validation_signed_no_top_nt():
    """signed=True with a meta that lacks 'top_nt' raises ValueError."""
    elist = pd.DataFrame(
        {"pre": [1, 2], "post": [2, 3], "count": [5, 6]}
    )
    meta_no_top_nt = pd.DataFrame({"root_id": [1, 2, 3]})
    with pytest.raises(ValueError, match="top_nt"):
        InfluenceCalculator.from_dataframes(
            elist, meta_df=meta_no_top_nt, signed=True,
            inhibitory_nts={'gaba'})


# ---------------------------------------------------------------------------
# 5b. excluded_nts removes outgoing edges and requires top_nt
# ---------------------------------------------------------------------------

def test_excluded_nts_removes_edges(edgelist_df, meta_df):
    """Pre-neurons whose top_nt is in excluded_nts contribute no edges
    to W. Comparing the dense column sums of W between the unfiltered
    and filtered builds shows zero outgoing weight for excluded neurons
    and unchanged outgoing weight for the rest."""
    pytest.importorskip("petsc4py", reason="PETSc not installed")

    excluded = {'glutamate', 'dopamine'}
    ic_all = InfluenceCalculator.from_dataframes(
        edgelist_df.copy(), meta_df, count_thresh=0)
    ic_excl = InfluenceCalculator.from_dataframes(
        edgelist_df.copy(), meta_df, count_thresh=0,
        excluded_nts=excluded)

    excluded_ids = set(meta_df.loc[meta_df['top_nt'].isin(excluded),
                                   'root_id'])
    excluded_idx = {ic_all.id_to_index[i] for i in excluded_ids
                    if i in ic_all.id_to_index}

    n = ic_all.n_neurons
    col_sum_all = np.zeros(n)
    col_sum_excl = np.zeros(n)
    for i in range(n):
        cols, vs = ic_all.W.getRow(i)
        for c, v in zip(cols, vs):
            col_sum_all[c] += v
        cols, vs = ic_excl.W.getRow(i)
        for c, v in zip(cols, vs):
            col_sum_excl[c] += v

    for j in range(n):
        if j in excluded_idx:
            assert abs(col_sum_excl[j]) < 1e-12, (
                f"excluded neuron col {j} should sum to 0")
        else:
            assert abs(col_sum_all[j] - col_sum_excl[j]) < 1e-9, (
                f"non-excluded col {j} changed: "
                f"{col_sum_all[j]} vs {col_sum_excl[j]}")


def test_excluded_nts_requires_top_nt():
    """excluded_nts without a top_nt column in meta raises ValueError."""
    elist = pd.DataFrame(
        {"pre": [1, 2], "post": [2, 3], "count": [5, 6]}
    )
    meta_no_top_nt = pd.DataFrame({"root_id": [1, 2, 3]})
    with pytest.raises(ValueError, match="top_nt"):
        InfluenceCalculator.from_dataframes(
            elist, meta_df=meta_no_top_nt,
            excluded_nts={'glutamate'})


# ---------------------------------------------------------------------------
# 6. Norm auto-computation
# ---------------------------------------------------------------------------

def test_norm_auto_computation(edgelist_df):
    """When 'norm' is absent, it is computed so each post neuron's weights
    sum to 1.0."""
    elist_no_norm = edgelist_df.drop(columns=["norm"])
    elist_out, _ = _validate_and_prepare_edgelist(
        elist_no_norm, None, False, 5
    )

    assert "norm" in elist_out.columns
    norm_sums = elist_out.groupby("post")["norm"].sum()
    assert np.allclose(norm_sums, 1.0), (
        "norm values do not sum to 1.0 per post neuron"
    )


# ---------------------------------------------------------------------------
# 7. Round-trip smoke test (requires PETSc + SLEPc)
# ---------------------------------------------------------------------------

def test_round_trip_smoke(edgelist_path, meta_path, meta_df):
    """Full pipeline: load → calculate_influence → adjust_influence."""
    pytest.importorskip("petsc4py", reason="PETSc not installed")
    pytest.importorskip("slepc4py", reason="SLEPc not installed")

    ic = InfluenceCalculator.from_csv(str(edgelist_path), str(meta_path))

    # Pick the first two neurons as seeds
    seed_ids = list(ic.id_to_index.keys())[:2]
    influence_df = ic.calculate_influence(seed_ids)

    assert len(influence_df) > 0
    score_col = "Influence_score_(unsigned)"
    assert score_col in influence_df.columns
    assert np.isfinite(influence_df[score_col]).all()

    result = adjust_influence(influence_df)
    assert len(result) > 0
    assert np.isfinite(result["adjusted_influence"]).all()
    assert (result["adjusted_influence"] >= 0).all()
