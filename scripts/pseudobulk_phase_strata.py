"""Cell-cycle-aware pseudobulks for the correlation analysis.

Two complementary views:

  strata  - one pseudobulk per (cell line, guide, phase). Correlations computed
            inside a stratum cannot be driven by phase composition, because
            every sample in it is the same phase.
  matched - every guide subsampled to a common phase composition within its
            cell line, so composition is equalised while all phases are kept.

A549 has only 4-14 G2M cells per guide, so its non-G1 strata are dropped by the
cell-count threshold; the matched view still covers both lines.
"""

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

H5AD = "/home/ubuntu/projects/symphony_10x_collateral/data/h5ad/SYM10X0_processed_filtered.h5ad"
OUT = "/home/ubuntu/projects/symphony_10x_collateral/output"
BLOCK = 2000
MIN_CELLS = 100
SEED = 1


def read_categorical(f, key):
    g = f[key]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
    codes = g["codes"][:]
    return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")


def sum_by_group(f, group_index, n_groups, n_genes, indptr, n_obs):
    """Sum raw counts per group, streaming the CSR matrix in row blocks."""
    sums = np.zeros((n_groups, n_genes), dtype=np.float64)
    counts = np.zeros(n_groups, dtype=np.int64)
    data = f["layers/counts/data"]
    indices = f["layers/counts/indices"]

    for start in range(0, n_obs, BLOCK):
        stop = min(start + BLOCK, n_obs)
        lo, hi = indptr[start], indptr[stop]
        block = sp.csr_matrix(
            (data[lo:hi], indices[lo:hi], indptr[start : stop + 1] - lo),
            shape=(stop - start, n_genes),
        )
        rows = np.arange(start, stop)
        sel = group_index[rows] >= 0
        if not sel.any():
            continue
        block = block[sel]
        groups = group_index[rows[sel]]
        for g in np.unique(groups):
            sums[g] += np.asarray(block[groups == g].sum(axis=0)).ravel()
            counts[g] += int((groups == g).sum())
    return sums, counts


f = h5py.File(H5AD, "r")

genes = np.array([g.decode() if isinstance(g, bytes) else g for g in f["var/_index"][:]])
guide = read_categorical(f, "obs/Hashsolo_classification")
cell_line = read_categorical(f, "obs/cell_line")
phase = read_categorical(f, "obs/phase")
indptr = f["layers/counts/indptr"][:]
n_obs = indptr.size - 1
n_genes = genes.size

assigned = ~np.isin(guide, ["Doublet", "Negative"])
expressed = pd.read_csv(f"{OUT}/pseudobulk_counts.csv", usecols=[0], index_col=0).index
gene_in_set = np.isin(genes, expressed)


def write_pseudobulk(labels, mask, path, metadata_path=None):
    levels = sorted(pd.unique(labels[mask]))
    lookup = {name: i for i, name in enumerate(levels)}
    index = np.full(n_obs, -1, dtype=np.int64)
    index[mask] = [lookup[name] for name in labels[mask]]

    sums, counts = sum_by_group(f, index, len(levels), n_genes, indptr, n_obs)
    table = pd.DataFrame(sums, index=levels, columns=genes)
    keep = counts >= MIN_CELLS
    print(f"{path}: {keep.sum()} of {len(levels)} pseudobulks pass {MIN_CELLS} cells")

    table.loc[keep, gene_in_set].T.astype(int).to_csv(f"{OUT}/{path}")
    if metadata_path is not None:
        pd.DataFrame(
            {"sample_id": np.array(levels)[keep], "n_cells": counts[keep]}
        ).to_csv(f"{OUT}/{metadata_path}", index=False)
    return pd.Series(counts, index=levels)


# --- phase strata -----------------------------------------------------------
strata_labels = np.array(
    [f"{c}|{g}|{p}" for c, g, p in zip(cell_line, guide, phase)], dtype=object
)
strata_counts = write_pseudobulk(
    strata_labels, assigned, "pseudobulk_counts_by_phase.csv",
    "pseudobulk_metadata_by_phase.csv",
)

# --- composition matching ---------------------------------------------------
# Target mix is the mRuby controls' average within each line; each guide is then
# subsampled to the largest cell count that reproduces that mix exactly.
rng = np.random.default_rng(SEED)
cells = pd.DataFrame(
    {"cell_line": cell_line, "guide": guide, "phase": phase, "row": np.arange(n_obs)}
)[assigned]

phases = sorted(cells["phase"].unique())
matched_rows = []

for line, line_cells in cells.groupby("cell_line"):
    controls = line_cells[line_cells["guide"].str.lower().str.startswith("mruby")]
    target = (
        controls["phase"].value_counts(normalize=True).reindex(phases).fillna(0)
    )
    print(f"\n{line} target phase mix: {target.round(3).to_dict()}")

    for guide_name, guide_cells in line_cells.groupby("guide"):
        available = guide_cells["phase"].value_counts().reindex(phases).fillna(0)
        usable = target[target > 0]
        total = int(np.floor(min(available[p] / target[p] for p in usable.index)))

        for p in usable.index:
            take = int(round(total * target[p]))
            pool = guide_cells.loc[guide_cells["phase"] == p, "row"].to_numpy()
            if take > 0:
                matched_rows.append(rng.choice(pool, size=min(take, pool.size), replace=False))
        print(f"  {guide_name}: {total} cells of {len(guide_cells)}")

matched = np.concatenate(matched_rows)
matched_mask = np.zeros(n_obs, dtype=bool)
matched_mask[matched] = True

matched_labels = np.array(
    [f"{c}|{g}" for c, g in zip(cell_line, guide)], dtype=object
)
write_pseudobulk(
    matched_labels, matched_mask, "pseudobulk_counts_phase_matched.csv",
    "pseudobulk_metadata_phase_matched.csv",
)

f.close()
