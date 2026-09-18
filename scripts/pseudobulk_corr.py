"""Pseudobulk correlation heatmaps for SYM10X0 Cas13d collateral-activity screen.

Cells are demultiplexed to a single guide via the HashSolo call, aggregated to
guide x cell_line pseudobulks from the raw counts layer, CPM+log1p normalised,
and correlated.  Two views are produced:
  1. expression correlation between pseudobulks
  2. correlation of per-guide logFC vs the mruby non-targeting controls, which
     is where a shared Cas13d collateral response would show up
"""

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns

H5AD = "/lustre/scratch125/casm/teams/team113/secure-lustre/team113/sym10x_corr_tmp/SYM10X0_processed_filtered.h5ad"
OUT = "/home/ubuntu/projects/symphony_10x_collateral/output"

MIN_CELLS = 50        # drop pseudobulks built from too few cells
MIN_MEAN_CPM = 1.0    # gene detection filter applied on the pseudobulk matrix


def read_categorical(f, key):
    g = f[key]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
    codes = g["codes"][:]
    out = np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")
    return out


f = h5py.File(H5AD, "r")

guide = read_categorical(f, "obs/Hashsolo_classification")
cell_line = read_categorical(f, "obs/cell_line")
genes = np.array([g.decode() if isinstance(g, bytes) else g for g in f["var/_index"][:]])

# HashSolo doublets/negatives carry no unambiguous guide assignment
keep = ~np.isin(guide, ["Doublet", "Negative"])
print(f"cells: {keep.sum()} singlet-assigned of {keep.size}")

phase = read_categorical(f, "obs/phase")

groups = pd.Series([f"{c}|{g}" for c, g in zip(cell_line[keep], guide[keep])])
group_levels = sorted(groups.unique())
gidx = groups.map({g: i for i, g in enumerate(group_levels)}).to_numpy()

# sum raw counts per group, streaming the CSR matrix in row blocks
n_genes = genes.size
sums = np.zeros((len(group_levels), n_genes), dtype=np.float64)
n_cells = np.zeros(len(group_levels), dtype=np.int64)

data = f["layers/counts/data"]
indices = f["layers/counts/indices"]
indptr = f["layers/counts/indptr"][:]
n_obs = indptr.size - 1

keep_idx = np.where(keep)[0]
pos_in_keep = np.full(n_obs, -1, dtype=np.int64)
pos_in_keep[keep_idx] = np.arange(keep_idx.size)

BLOCK = 2000
for start in range(0, n_obs, BLOCK):
    stop = min(start + BLOCK, n_obs)
    lo, hi = indptr[start], indptr[stop]
    block = sp.csr_matrix(
        (data[lo:hi], indices[lo:hi], indptr[start:stop + 1] - lo),
        shape=(stop - start, n_genes),
    )
    rows = np.arange(start, stop)
    sel = pos_in_keep[rows] >= 0
    if not sel.any():
        continue
    block = block[sel]
    g_block = gidx[pos_in_keep[rows[sel]]]
    for g in np.unique(g_block):
        sums[g] += np.asarray(block[g_block == g].sum(axis=0)).ravel()
        n_cells[g] += int((g_block == g).sum())

f.close()

pb = pd.DataFrame(sums, index=group_levels, columns=genes)
pb["__n_cells"] = n_cells
pb = pb[pb["__n_cells"] >= MIN_CELLS]
n_cells_kept = pb.pop("__n_cells")
print(pd.Series(n_cells_kept.values, index=pb.index, name="cells").to_string())

# CPM + log1p
cpm = pb.div(pb.sum(axis=1), axis=0) * 1e6
expressed = cpm.mean(axis=0) >= MIN_MEAN_CPM
logcpm = np.log1p(cpm.loc[:, expressed])
print(f"genes retained: {expressed.sum()} of {expressed.size}")
logcpm.T.to_csv(f"{OUT}/pseudobulk_logcpm.csv")

# raw summed counts and cell numbers, for count-based DE downstream
pb.loc[:, expressed].T.astype(int).to_csv(f"{OUT}/pseudobulk_counts.csv")
# cell cycle composition per pseudobulk: targeting guides shift the phase mix,
# so downstream models need it as a covariate rather than reading it as signal
phase_fractions = (
    pd.crosstab(groups.to_numpy(), phase[keep], normalize="index")
    .reindex(pb.index)
    .add_prefix("frac_")
)

pd.DataFrame({
    "sample_id": pb.index,
    "cell_line": [i.split("|")[0] for i in pb.index],
    "guide": [i.split("|")[1] for i in pb.index],
    "n_cells": n_cells_kept.values,
}).join(phase_fractions, on="sample_id").to_csv(
    f"{OUT}/pseudobulk_metadata.csv", index=False
)

# G1-only pseudobulks, as a cell-cycle-matched sensitivity analysis
g1 = keep & (phase == "G1")
g1_groups = pd.Series([f"{c}|{g}" for c, g in zip(cell_line[g1], guide[g1])])
g1_levels = sorted(g1_groups.unique())
g1_gidx = g1_groups.map({g: i for i, g in enumerate(g1_levels)}).to_numpy()
g1_sums = np.zeros((len(g1_levels), n_genes), dtype=np.float64)
g1_n = np.zeros(len(g1_levels), dtype=np.int64)

g1_idx = np.where(g1)[0]
pos_in_g1 = np.full(n_obs, -1, dtype=np.int64)
pos_in_g1[g1_idx] = np.arange(g1_idx.size)

with h5py.File(H5AD, "r") as fh:
    data_g1 = fh["layers/counts/data"]
    indices_g1 = fh["layers/counts/indices"]
    for start in range(0, n_obs, BLOCK):
        stop = min(start + BLOCK, n_obs)
        lo, hi = indptr[start], indptr[stop]
        block = sp.csr_matrix(
            (data_g1[lo:hi], indices_g1[lo:hi], indptr[start:stop + 1] - lo),
            shape=(stop - start, n_genes),
        )
        rows = np.arange(start, stop)
        sel = pos_in_g1[rows] >= 0
        if not sel.any():
            continue
        block = block[sel]
        g_block = g1_gidx[pos_in_g1[rows[sel]]]
        for g in np.unique(g_block):
            g1_sums[g] += np.asarray(block[g_block == g].sum(axis=0)).ravel()
            g1_n[g] += int((g_block == g).sum())

g1_pb = pd.DataFrame(g1_sums, index=g1_levels, columns=genes)
g1_pb = g1_pb[g1_n >= MIN_CELLS]
g1_pb.loc[:, expressed].T.astype(int).to_csv(f"{OUT}/pseudobulk_counts_g1.csv")
print(f"G1-only pseudobulks: {g1_pb.shape[0]} passing {MIN_CELLS} cells")


def heatmap(corr, title, fname, cmap, center, row_colors=None):
    g = sns.clustermap(
        corr, cmap=cmap, center=center, figsize=(9, 8.5),
        annot=True, fmt=".2f", annot_kws={"size": 6},
        row_colors=row_colors, col_colors=row_colors,
        cbar_kws={"label": "Pearson r"}, linewidths=0.4,
    )
    g.ax_heatmap.set_xlabel("")
    g.ax_heatmap.set_ylabel("")
    g.figure.suptitle(title, y=1.0, fontsize=11)
    g.savefig(fname, dpi=200, bbox_inches="tight")
    plt.close(g.figure)


lines = pd.Series([i.split("|")[0] for i in logcpm.index], index=logcpm.index)
targets = pd.Series([i.split("|")[1].split("_")[0].upper() for i in logcpm.index], index=logcpm.index)
line_colors = lines.map(dict(zip(lines.unique(), sns.color_palette("Set2", lines.nunique()))))

corr_expr = logcpm.T.corr(method="pearson")
corr_expr.to_csv(f"{OUT}/correlation_expression.csv")
heatmap(corr_expr, "Pseudobulk expression correlation (log1p CPM)",
        f"{OUT}/heatmap_expression_correlation.png", "viridis", None, line_colors)

# logFC vs the mruby non-targeting controls, computed within each cell line.
# A control is compared against the *other* controls only: including a sample in
# its own baseline forces the controls to correlate at -1/(n-1) by construction.
fc_cols = {}
for line in lines.unique():
    members = logcpm.index[lines == line]
    ctrl = [m for m in members if m.split("|")[1].lower().startswith("mruby")]
    if not ctrl:
        print(f"no mruby control for {line}, skipping")
        continue
    for m in members:
        ref = [c for c in ctrl if c != m]
        fc_cols[m] = logcpm.loc[m] - logcpm.loc[ref].mean(axis=0)

fc = pd.DataFrame(fc_cols).T
fc.T.to_csv(f"{OUT}/pseudobulk_logfc_vs_mruby.csv")
corr_fc = fc.T.corr(method="pearson")
corr_fc.to_csv(f"{OUT}/correlation_logfc.csv")
heatmap(corr_fc, "Correlation of logFC vs mruby controls\n(shared transcriptional response to Cas13d activity)",
        f"{OUT}/heatmap_logfc_correlation.png", "RdBu_r", 0, line_colors)

print("\nwritten to", OUT)


# Consensus collateral signature: genes moving the same way across all targeting
# guides within a cell line, scored as mean logFC / SD across guides.
targeting = [i for i in fc.index if not i.split("|")[1].lower().startswith("mruby")]
sig = {}
for line in lines.unique():
    members = [i for i in targeting if i.startswith(line + "|")]
    if len(members) < 3:
        continue
    sub = fc.loc[members]
    sig[f"{line}_mean_logFC"] = sub.mean(axis=0)
    sig[f"{line}_consistency"] = sub.mean(axis=0) / sub.std(axis=0)
sig = pd.DataFrame(sig).replace([np.inf, -np.inf], np.nan).dropna()
# require a real effect in both lines, else consistency scores tiny noise highly
mean_cols = [c for c in sig if c.endswith("mean_logFC")]
sig = sig[(sig[mean_cols].abs() >= 0.05).all(axis=1)]
sig["rank_score"] = sig[[c for c in sig if c.endswith("consistency")]].mean(axis=1)
sig.sort_values("rank_score").to_csv(f"{OUT}/consensus_collateral_signature.csv")
print(f"\nsignature genes passing |logFC| >= 0.05 in both lines: {len(sig)}")
print("\ntop consistently DOWN across targeting guides:")
print(sig.sort_values("rank_score").head(15).round(3).to_string())
print("\ntop consistently UP across targeting guides:")
print(sig.sort_values("rank_score").tail(15).round(3).to_string())
