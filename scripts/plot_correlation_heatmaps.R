library(here)
source(here("scripts", "heatmap_helpers.R"))

# Pseudobulk correlation heatmaps for the SYM10X0 Cas13d collateral screen.
# Matrices come from scripts/pseudobulk_corr.py; sample ids are "cell_line|guide".

corr_logfc <- read_corr(here("output", "correlation_logfc.csv"))
corr_expression <- read_corr(here("output", "correlation_expression.csv"))

# mRuby guides define the logFC baseline, so their correlations there are fixed
# by construction (controls at -1/(n-1), controls vs targeting at ~0) and carry
# no signal; they stay in the expression panel, where they are ordinary samples
targeting_samples <- rownames(corr_logfc) |>
  str_subset("mruby", negate = TRUE) |>
  discard(is_ineffective)
corr_logfc <- corr_logfc[targeting_samples, targeting_samples]

annotation_logfc <- build_annotation_table(rownames(corr_logfc))
annotation_expression <- build_annotation_table(rownames(corr_expression))

# palettes are built from the full sample set so colours mean the same thing in
# both panels, even though mRuby is only shown in the expression one
annotation_colours <- list(
  `Cell line` = make_palette(annotation_expression[["cell_line"]], "Dark2"),
  Target = make_palette(annotation_expression[["target"]], "Set2"),
  Guide = make_palette(annotation_expression[["guide"]], "Paired")
)

# logFC correlations straddle zero, so use a diverging scale centred on white;
# expression correlations are all high and positive, so use a sequential one
logfc_limit <- max(abs(corr_logfc[upper.tri(corr_logfc)]))
logfc_colours <- colorRamp2(
  c(-logfc_limit, 0, logfc_limit),
  c("#2166AC", "white", "#B2182B")
)
expression_colours <- colorRamp2(
  range(corr_expression[upper.tri(corr_expression)]),
  c("#2D004B", "#FDE725")
)

plot_correlation_heatmap(
  corr_logfc,
  annotation_logfc,
  "Correlation of logFC vs mRuby controls\nshared transcriptional response to Cas13d activity",
  logfc_colours,
  "Pearson r",
  annotation_colours,
  legend_breaks = seq(-round(logfc_limit, 2), round(logfc_limit, 2), length.out = 5) |> round(2)
) |>
  save_heatmap(here("output", "heatmap_logfc_correlation_complex.png"))

plot_correlation_heatmap(
  corr_expression,
  annotation_expression,
  "Pseudobulk expression correlation (log1p CPM)",
  expression_colours,
  "Pearson r",
  annotation_colours,
  legend_breaks = round(seq(min(corr_expression[upper.tri(corr_expression)]), max(corr_expression[upper.tri(corr_expression)]), length.out = 5), 2)
) |>
  save_heatmap(here("output", "heatmap_expression_correlation_complex.png"))

save_per_line_panels(
  corr_logfc, "heatmap_logfc_correlation", "logFC vs mRuby controls",
  diverging = TRUE, annotation_colours = annotation_colours
)
save_per_line_panels(
  corr_expression, "heatmap_expression_correlation", "pseudobulk expression",
  diverging = FALSE, annotation_colours = annotation_colours
)
