"""CD274 counts per G1 pseudobulk.

CD274 averages ~0.8 CPM and falls below the expression filter applied to the
pseudobulk matrices, so its counts are read back from the h5ad directly. The
raw totals are kept in the output so downstream code can flag how few counts
the estimate rests on (1-9 per pseudobulk).
"""

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
column = int(np.where(genes == "CD274")[0][0])

guide = read_categorical(f, "obs/Hashsolo_classification")
cell_line = read_categorical(f, "obs/cell_line")
phase = read_categorical(f, "obs/phase")
keep = (~np.isin(guide, ["Doublet", "Negative"])) & (phase == "G1")

indptr = f["layers/counts/indptr"][:]
n_obs = indptr.size - 1
cd274 = np.zeros(n_obs)
total = np.zeros(n_obs)
data = f["layers/counts/data"]
indices = f["layers/counts/indices"]

for start in range(0, n_obs, BLOCK):
    stop = min(start + BLOCK, n_obs)
    lo, hi = indptr[start], indptr[stop]
    block = sp.csr_matrix(
        (data[lo:hi], indices[lo:hi], indptr[start : stop + 1] - lo),
        shape=(stop - start, genes.size),
    )
    cd274[start:stop] = np.asarray(block[:, column].todense()).ravel()
    total[start:stop] = np.asarray(block.sum(axis=1)).ravel()

cells = pd.DataFrame(
    {"cell_line": cell_line, "guide": guide, "cd274_counts": cd274, "total_counts": total}
)[keep]

summary = (
    cells.groupby(["cell_line", "guide"])
    .agg(
        cd274_counts=("cd274_counts", "sum"),
        total_counts=("total_counts", "sum"),
        n_cells=("cd274_counts", "size"),
        n_detected=("cd274_counts", lambda x: (x > 0).sum()),
    )
    .reset_index()
)
summary["cd274_cpm"] = 1e6 * summary["cd274_counts"] / summary["total_counts"]
summary.to_csv(f"{OUT}/cd274_counts_g1.csv", index=False)
print(summary.to_string(index=False))
f.close()
