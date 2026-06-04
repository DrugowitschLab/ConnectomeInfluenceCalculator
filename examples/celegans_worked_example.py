"""Worked example: influence of all sensory neurons on every other
neuron in the C. elegans connectome.

Loads the bundled C. elegans dataset, computes per-seed influence
from every sensory neuron onto every non-sensory target, collapses
bilateral pairs into cell classes by stripping the trailing L/R
suffix on both seed and target sides (so AVAL/AVAR → AVA, AVDL/AVDR →
AVD, IL2DL/IL2DR → IL2D), sums the raw influence per
(target_class, seed_class), log-adjusts via adjust_influence with an
auto-calibrated const, and renders the result as two heatmaps
(unsigned and signed).

The heatmaps show the raw adjusted_influence values directly — no
per-row min-max rescaling is applied — anchored at 0 and at the data
extremum of each variant. Rows (seed classes) and columns (target
classes) are grouped by anatomical body_part and ordered by
average-linkage hierarchical clustering within each group.

Run from the repository root:

    python examples/celegans_worked_example.py
"""
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.cluster.hierarchy import linkage, leaves_list
from tqdm import tqdm

from InfluenceCalculator import InfluenceCalculator
from InfluenceCalculator.data import celegans_edgelist, celegans_meta


REPO_ROOT = Path(__file__).resolve().parent.parent
OUT_DIR = REPO_ROOT / "docs" / "images"

# Neurotransmitter handling for C. elegans (signed mode only).
#   * Acetylcholine — dominant excitatory transmitter; left positive.
#   * Glutamate     — net-excitatory at most postsynaptic partners
#                     (AMPA-/NMDA-like GLR-* / NMR-* receptors), even
#                     though GluCl chloride channels make it inhibitory
#                     at a minority of synapses (e.g. AWC → AIY). We
#                     leave it positive by default; recovering the
#                     minority inhibitory cases needs edge-level
#                     overrides, which the NT-level API does not
#                     currently support.
#   * GABA          — dominant inhibitory transmitter; negated below
#                     via inhibitory_nts.
#   * Dopamine, serotonin, octopamine — neuromodulators whose sign at
#                     a given target depends on the receptor mix;
#                     silence their pre-neurons via excluded_nts.
#
# The library has no per-organism defaults; both sets are user input.
# Two reasonable starting points:
#
#   Drosophila (the BANC pipeline's historical convention):
#     inhibitory_nts = {'glutamate', 'gaba', 'serotonin', 'octopamine'}
#     excluded_nts   = set()
#
#   C. elegans (only ACh and GABA have unambiguous signs):
#     inhibitory_nts = {'gaba'}
#     excluded_nts   = {'glutamate', 'dopamine', 'serotonin', 'octopamine'}
#
# Note that the C. elegans signed configuration excludes the majority
# of sensory drivers: of the 83 non-pharyngeal sensory neurons, 47 are
# glutamatergic, 8 dopaminergic, 2 serotonergic, 13 cholinergic, and 13
# unannotated — only 26 actually transmit in the signed heatmap, and
# the silenced classes appear as a wide blank band on the seed axis.
# The unsigned variant has no such exclusion: every sensory neuron
# drives a column.
CELEGANS_INHIBITORY_NTS = {'gaba'}
CELEGANS_EXCLUDED_NTS = {'dopamine', 'serotonin', 'octopamine'}

# Seeds are sensory neurons; targets are everything else (interneurons,
# motor neurons, and modulatory neurons). Pharyngeal neurons are
# excluded because the pharyngeal nervous system is essentially
# isolated from the rest of the connectome and would otherwise
# dominate the heatmap with a block of zeros.
SEED_CLASS = 'sensory'
EXCLUDE_BODY_PARTS = {'pharynx'}

# Minimum synaptic count for an edge to be retained (i.e. require count
# >= COUNT_THRESH). With count_thresh=0 every weak / single-synapse edge
# is retained.
COUNT_THRESH = 0

# Target spectral radius of the rescaled W̃.  The amplification of the
# leading eigenmode in (I - W̃)^-1 is 1 / (1 - lambda_max), so the
# library default 0.99 (~100x) makes the leading mode dominate and
# every column of (I - W̃)^-1 share the same shape.  Setting it to 0.5
# (~2x) damps that mode and exposes per-target seed-specificity, at
# the cost of attenuating long polysynaptic paths.
LAMBDA_MAX = 0.5


def build_calculator(edges, meta, signed):
    """Construct an InfluenceCalculator at the configured threshold.

    Unsigned mode is a structural propagation view: every edge counts
    positively, regardless of its presynaptic neurotransmitter, so
    neither inhibitory_nts nor excluded_nts is applied.

    Signed mode applies the biology: GABA pre-neurons are negated, and
    pre-neurons whose top_nt cannot be assigned a single sign safely
    (glutamate, dopamine, serotonin, octopamine — gating both excitatory
    and inhibitory receptors in C. elegans) are excluded entirely.
    """
    return InfluenceCalculator(
        edges, meta, signed=signed, count_thresh=COUNT_THRESH,
        inhibitory_nts=(CELEGANS_INHIBITORY_NTS if signed else None),
        excluded_nts=(CELEGANS_EXCLUDED_NTS if signed else None),
        lambda_max=LAMBDA_MAX)


def per_seed_influence(ic, seed_ids, score_col, desc):
    """Run calculate_influence once per seed and concatenate the results
    into a long-form DataFrame with columns target, seed, is_seed and
    the raw influence score.
    """
    rows = []
    for seed_id in tqdm(seed_ids, desc=desc):
        # Skip the per-row adjusted columns -- they are recomputed
        # below at the (target_class, seed_class) granularity.
        df = ic.calculate_influence([seed_id], adjust=False)
        df = df.rename(columns={'id': 'target'})
        df['seed'] = seed_id
        rows.append(df[['target', 'seed', 'is_seed', score_col]])
    return pd.concat(rows, ignore_index=True)


def cluster_within_groups(matrix, group_of):
    """Return (ordering, boundaries) where ordering is a list of
    matrix.index labels grouped by group_of (visited alphabetically)
    and ordered by average-linkage hierarchical clustering within each
    group, and boundaries is a list of (group_label, cumulative_index)
    tuples marking the right edge of each group.
    """
    groups = pd.Series(matrix.index).map(group_of)
    order = []
    boundaries = []
    cumulative = 0
    for group in sorted(groups.dropna().unique()):
        members = matrix.index[groups.values == group]
        if len(members) == 1:
            order.extend(members)
        else:
            sub = matrix.loc[members].fillna(0).values
            link = linkage(sub, method='average')
            order.extend(np.array(members)[leaves_list(link)])
        cumulative += len(members)
        boundaries.append((group, cumulative))
    return order, boundaries


def cell_class(name):
    """Strip the trailing L/R from a C. elegans neuron name to collapse
    bilaterally / quadripartite-symmetric pairs into one cell class
    (AVAL/AVAR -> AVA, AVDL/AVDR -> AVD, IL2DL/IL2DR -> IL2D).
    Numerically suffixed motor neurons (VD01, VD13, ...) and singletons
    are returned unchanged.

    The lookbehind requires an upper-case letter so VD01 / AS11 are not
    truncated.  We deliberately do NOT include DL/DR/VL/VR as
    alternatives: leftmost matching would incorrectly turn AVDL into
    AV (matching DL at position 2) instead of AVD.  Stripping just the
    final L|R gives the right answer for both AVDL→AVD and IL2DL→IL2D.
    """
    return re.sub(r'(?<=[A-Z])[LR]$', '', name)


def render_heatmap(matrix, row_bp, col_bp, title, out_path, signed,
                   row_label, col_label):
    """Render a (row x col) heatmap of adjusted_influence with rows and
    columns grouped by body_part and clustered within each group. The
    unsigned variant uses a sequential greyscale (low to high in
    [0, max]) and the signed variant uses a diverging blue→red
    (RdBu_r, white at 0) over [-bound, +bound] where bound = max |x|,
    so the two remain visually distinct and net inhibition / net
    excitation are symmetric in the signed view.

    row_bp and col_bp are pd.Series mapping root_id → body_part for the
    row and column indices respectively.
    """
    row_order, r_bounds = cluster_within_groups(matrix, row_bp)
    col_order, c_bounds = cluster_within_groups(matrix.T, col_bp)
    M = matrix.loc[row_order, col_order]

    if signed:
        bound = float(np.nanmax(np.abs(M.values)))
        cmap, vmin, vmax = 'RdBu_r', -bound, bound
    else:
        cmap = 'Greys'
        vmin, vmax = 0.0, float(np.nanmax(M.values))

    n_rows, n_cols = len(row_order), len(col_order)
    fig_w = max(11, 0.13 * n_cols + 4)
    fig_h = max(7, 0.16 * n_rows + 3)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    sns.heatmap(M, ax=ax, cmap=cmap, vmin=vmin, vmax=vmax,
                cbar_kws={'label': 'adjusted_influence'},
                xticklabels=True, yticklabels=True, linewidths=0)

    ax.set_title(title, fontsize=12)
    ax.set_xlabel(col_label, labelpad=22)
    ax.set_ylabel(row_label, labelpad=28)
    ax.tick_params(axis='x', labelsize=5, rotation=90)
    ax.tick_params(axis='y', labelsize=6)

    # Body-part group separators and labels along each axis. Group
    # labels are placed in axes coordinates a fixed distance outside the
    # heatmap so they clear the rotated tick labels regardless of how
    # long those labels are.
    prev = 0
    for group, cum in r_bounds:
        if cum < n_rows:
            ax.axhline(cum, color='orange', lw=1.0)
        y_data = (prev + cum) / 2
        y_axes = ax.transAxes.inverted().transform(
            ax.transData.transform((0, y_data)))[1]
        ax.text(-0.06, y_axes, group, transform=ax.transAxes,
                ha='right', va='center', fontsize=9,
                fontweight='bold', clip_on=False)
        prev = cum
    prev = 0
    for group, cum in c_bounds:
        if cum < n_cols:
            ax.axvline(cum, color='orange', lw=1.0)
        x_data = (prev + cum) / 2
        x_axes = ax.transAxes.inverted().transform(
            ax.transData.transform((x_data, 0)))[0]
        ax.text(x_axes, -0.12, group, transform=ax.transAxes,
                ha='center', va='top', fontsize=9,
                fontweight='bold', clip_on=False)
        prev = cum

    fig.tight_layout()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=130, bbox_inches='tight')
    plt.close(fig)
    print(f"wrote {out_path.relative_to(REPO_ROOT)}")


def main():
    edges = celegans_edgelist()
    meta = celegans_meta()

    not_excluded = ~meta['body_part'].isin(EXCLUDE_BODY_PARTS)
    seed_meta = meta[(meta['super_class'] == SEED_CLASS) & not_excluded]
    target_meta = meta[(meta['super_class'] != SEED_CLASS) & not_excluded]
    seed_ids = seed_meta['root_id'].tolist()
    target_ids = set(target_meta['root_id'])

    # Map every neuron to its cell class (AVAL/AVAR -> AVA), and pick
    # one body_part per class (the first member's; in practice all
    # members of a class share a body_part).
    name_to_class = {n: cell_class(n) for n in meta['root_id']}
    class_bp = (meta.assign(cell_class=meta['root_id'].map(name_to_class))
                    .drop_duplicates('cell_class')
                    .set_index('cell_class')['body_part'])

    print(f"{len(seed_ids)} sensory seeds ({seed_meta['root_id']
                                            .map(name_to_class)
                                            .nunique()} classes) → "
          f"{len(target_ids)} non-sensory targets "
          f"({target_meta['root_id'].map(name_to_class).nunique()} classes), "
          f"excluding {sorted(EXCLUDE_BODY_PARTS)}")

    for signed, label in [(False, 'unsigned'), (True, 'signed')]:
        ic = build_calculator(edges, meta, signed=signed)
        score_col = f'Influence_score_({label})'
        long_df = per_seed_influence(ic, seed_ids, score_col,
                                     desc=f'{label} influence')

        # Collapse bilateral pairs by replacing target / seed root_ids
        # with their cell class on both axes.  adjust_influence sums
        # the raw influence per (target_class, seed_class) group.
        sub_raw = long_df[long_df['target'].isin(target_ids)].copy()
        sub_raw['target'] = sub_raw['target'].map(name_to_class)
        sub_raw['seed'] = sub_raw['seed'].map(name_to_class)

        # Auto-calibrate const from the cell-class summed magnitudes so
        # the smallest non-zero |sum| just clips to 0.
        per_pair = (sub_raw.groupby(['target', 'seed'])[score_col]
                           .sum().abs())
        per_pair = per_pair[per_pair > 0]
        const = float(-np.log(per_pair.min()))

        adjusted = InfluenceCalculator.adjust_influence(sub_raw, const=const)
        matrix = adjusted.pivot_table(index='target', columns='seed',
                                      values='adjusted_influence')

        # Transpose so seed classes index the rows.
        matrix = matrix.T

        title = (f'C. elegans sensory → non-sensory influence — {label} '
                 f'(count_thresh={COUNT_THRESH}, lambda_max={LAMBDA_MAX}, '
                 f'const={const:.2f}, cell-class sum)')
        render_heatmap(
            matrix, class_bp, class_bp, title,
            OUT_DIR / f'influence_heatmap_{label}.png',
            signed=signed,
            row_label=f'sensory seed class (n = {matrix.shape[0]})',
            col_label=f'non-sensory target class '
                      f'(n = {matrix.shape[1]})')

        n_neg = int((adjusted['adjusted_influence'] < 0).sum())
        print(f"  {label}: {len(adjusted)} (target_class, seed_class) "
              f"cells, const={const:.2f}, "
              f"{n_neg} with negative adjusted_influence")


if __name__ == "__main__":
    main()
