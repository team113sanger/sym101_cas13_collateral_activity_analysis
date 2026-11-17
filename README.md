# SYM10X0 — Cas13d collateral activity

Analysis of a Cas13d perturb-seq screen in A549 and NCI-H1299, testing whether
guide-directed Cas13d produces a target-independent transcriptional change.

Input is `SYM10X0_processed_filtered.h5ad` (28,492 cells × 42,350 genes), already
demultiplexed by SNP donor and by guide (HashSolo) upstream. `data/h5ad` is a
symlink to the unpacked file.

## Running

R 4.4.0 via `rig`, dependencies pinned in `renv.lock`. Run in order; each step
reads the previous step's output from `output/`.

```
python3 scripts/pseudobulk_corr.py          # pseudobulks, phase fractions, correlations
python3 scripts/pseudobulk_phase_strata.py  # per-phase and phase-matched pseudobulks
python3 scripts/extract_cd274_counts.py     # CD274 sits below the expression filter

rig run -r 4.4.0 -f scripts/differential_expression.R   # DE + GSEA tables
rig run -r 4.4.0 -f scripts/correlation_by_phase.R      # FIGURE 1
rig run -r 4.4.0 -f scripts/figure_lfc_histogram.R      # FIGURE 2
rig run -r 4.4.0 -f scripts/figure_nes_concordance.R    # FIGURE 3

rig run -r 4.4.0 -f scripts/transcript_mass_budget.R    # numbers quoted in the text
rig run -r 4.4.0 -f scripts/rp_module_permutation.R     # module selection control
rig run -r 4.4.0 -f scripts/depmap_essentiality.R       # target essentiality (needs DepMap data)
```

## Figures

| Figure | File | Script |
|---|---|---|
| 1 — guide–guide correlation, G1 cells | `figure_guide_correlation_g1only_by_cell_line.pdf` | `correlation_by_phase.R` |
| 2 — per-gene logFC distribution, CD9_1 | `figure_target_vs_rp_histogram.pdf` | `figure_lfc_histogram.R` |
| 3 — pathway NES concordance between lines | `figure_nes_concordance.pdf` | `figure_nes_concordance.R` |

## Analysis decisions that affect every result

- **Fold changes are leave-one-out** against the mRuby controls of the same cell
  line; a control is never part of its own baseline, which would otherwise force
  the controls to correlate at −1/(n−1).
- **Cross-line comparisons use G1 cells only.** Targeting guides shift the cell
  cycle distribution in H1299 (G1 37.1% → 29.5%), and without this an E2F/G2M
  composition shift reads as an expression signature.
- **CD9_3 is excluded** from group-level statistics. It achieves −0.13 log2
  against a target measured at >130 CPM in both lines, so it is a non-functional
  guide. Defined once as `ineffective_guides` in `heatmap_helpers.R`.
- **Only CD9 has a trustworthy on-target readout.** FAS is 3.7 CPM in H1299 and
  CD274 ~0.8 CPM (1–9 counts per pseudobulk), so knockdown cannot be quantified
  for those targets in those lines.
- **The expression filter is re-applied per analysis**, so gene counts differ by
  cell subset: 15,983 across all cells, 15,788 in G1 only. The 195 genes lost are
  enriched ~8-fold for S/G2M-specific transcripts, which is expected.

## Limitations

The dataset cannot establish collateral RNase activity. There is no
catalytically-dead or enzyme-free arm, so the contrast is guide-that-matches
versus guide-that-does-not. Library-size normalisation measures relative
composition, so a uniform loss of transcripts is invisible and what is observed
is a redistribution. Effect sizes carry selection bias: the ribosomal module was
chosen after inspecting genome-wide enrichment results, with no held-out data.

## Authors

This analysis was performed by Jamie Billington (<jb63@sanger.ac.uk>) in the
Adams group at the Wellcome Sanger Institute.

## Contact

- Jamie Billington (<jb63@sanger.ac.uk>)

## Licence

This repository is released under the GNU Affero General Public License v3.0.
See [LICENSE](LICENSE) for the full text.
