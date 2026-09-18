"""Per-cell table for the single-cell dose analysis.

For every HashSolo-assigned cell: the ribosomal protein module score, the UMI
count of the guide it was assigned to, Cas13d expression, and the technical
covariates needed to control the regression.
"""

import re

import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

H5AD = "/home/ubuntu/projects/symphony_10x_collateral/data/h5ad/SYM10X0_processed_filtered.h5ad"
OUT = "/home/ubuntu/projects/symphony_10x_collateral/output"
BLOCK = 2000


def read_categorical(f, key):
    g = f[key]
    cats = np.array([c.decode() if isinstance(c, bytes) else c for c in g["categories"][:]])
    codes = g["codes"][:]
    return np.where(codes >= 0, cats[np.clip(codes, 0, None)], "NA")


f = h5py.File(H5AD, "r")

genes = np.array([g.decode() if isinstance(g, bytes) else g for g in f["var/_index"][:]])
guide = read_categorical(f, "obs/Hashsolo_classification")
cell_line = read_categorical(f, "obs/cell_line")
phase = read_categorical(f, "obs/phase")

# Same module as the pseudobulk analysis: the bare regex also catches RP
# pseudogenes, so intersect it with the expressed gene set used there.
expressed = set(
    pd.read_csv(f"{OUT}/pseudobulk_counts.csv", usecols=[0], index_col=0).index
)
rp_mask = np.array(
    [bool(re.match(r"^RP[LS]\d", g)) and g in expressed for g in genes]
)
mito_mask = np.array([g.startswith("MT-") for g in genes])
print(f"ribosomal protein genes: {rp_mask.sum()}, mitochondrial: {mito_mask.sum()}")

indptr = f["layers/counts/indptr"][:]
n_obs = indptr.size - 1
n_genes = genes.size

rp_counts = np.zeros(n_obs)
mito_counts = np.zeros(n_obs)
total_counts = np.zeros(n_obs)
n_detected = np.zeros(n_obs, dtype=np.int64)

data = f["layers/counts/data"]
indices = f["layers/counts/indices"]
for start in range(0, n_obs, BLOCK):
    stop = min(start + BLOCK, n_obs)
    lo, hi = indptr[start], indptr[stop]
    block = sp.csr_matrix(
        (data[lo:hi], indices[lo:hi], indptr[start : stop + 1] - lo),
        shape=(stop - start, n_genes),
    )
    rp_counts[start:stop] = np.asarray(block[:, rp_mask].sum(axis=1)).ravel()
    mito_counts[start:stop] = np.asarray(block[:, mito_mask].sum(axis=1)).ravel()
    total_counts[start:stop] = np.asarray(block.sum(axis=1)).ravel()
    n_detected[start:stop] = block.getnnz(axis=1)

# UMI count of the guide each cell was actually assigned to
guide_capture = pd.DataFrame(
    {
        col: f[f"obsm/guide_capture/{col}"][:]
        for col in f["obsm/guide_capture"]
        if col != "cellbarcode"
    }
)
assigned_umi = np.array(
    [
        guide_capture.at[i, g] if g in guide_capture.columns else np.nan
        for i, g in enumerate(guide)
    ],
    dtype=float,
)

targeting_cols = [c for c in guide_capture.columns if not c.lower().startswith("mruby")]
control_cols = [c for c in guide_capture.columns if c.lower().startswith("mruby")]

cells = pd.DataFrame(
    {
        "cell_line": cell_line,
        "guide": guide,
        "phase": phase,
        "cas13_expr": f["obs/Cas13_expr"][:],
        "n_unique_guides": f["obs/n_unique_guides"][:],
        "guide_umi_total": f["obs/guide_UMI_counts"][:],
        "assigned_guide_umi": assigned_umi,
        "targeting_guide_umi": guide_capture[targeting_cols].sum(axis=1).to_numpy(),
        "control_guide_umi": guide_capture[control_cols].sum(axis=1).to_numpy(),
        "rp_fraction": rp_counts / total_counts,
        "mito_fraction": mito_counts / total_counts,
        "total_counts": total_counts,
        "n_genes_detected": n_detected,
    }
)

cells = cells[~cells["guide"].isin(["Doublet", "Negative"])].copy()
cells["condition"] = np.where(
    cells["guide"].str.lower().str.startswith("mruby"), "non_targeting", "targeting"
)
cells["target"] = (
    cells["guide"].str.upper().str.replace(r"_\d+$", "", regex=True).replace({"CD95": "FAS"})
)

cells.to_csv(f"{OUT}/single_cell_dose.csv", index=False)
print(f"cells written: {len(cells)}")
print(cells.groupby(["cell_line", "condition"])[["rp_fraction", "assigned_guide_umi", "cas13_expr"]].mean().round(4))
f.close()
