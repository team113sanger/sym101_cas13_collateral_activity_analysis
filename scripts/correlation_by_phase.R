library(here)
source(here("scripts", "heatmap_helpers.R"))

# Cell-cycle-aware versions of the guide-guide correlations.
#
#   strata  - correlations computed within a single phase, so phase composition
#             cannot contribute at all
#   matched - all phases retained, but every guide subsampled to a common phase
#             mix within its cell line
#
# Both read pseudobulks from scripts/pseudobulk_phase_strata.py.

min_pairs <- 3

strata_counts <- read_corr(here("output", "pseudobulk_counts_by_phase.csv"))
matched_counts <- read_corr(here("output", "pseudobulk_counts_phase_matched.csv"))

# --- correlations within each phase -----------------------------------------
strata <- tibble(sample_id = colnames(strata_counts)) |>
  separate_wider_delim(
    sample_id,
    delim = "|", names = c("cell_line", "guide", "phase"), cols_remove = FALSE
  )

# a stratum needs its own mRuby controls to form logFCs, and enough targeting
# guides left over for the correlation to mean anything
usable <- strata |>
  summarise(
    n_controls = sum(str_detect(str_to_lower(guide), "mruby")),
    n_targeting = sum(!str_detect(str_to_lower(guide), "mruby")),
    .by = c(cell_line, phase)
  ) |>
  filter(n_controls >= 2, n_targeting >= min_pairs)

cat("=== usable phase strata ===\n")
print(usable)

phase_correlations <- usable |>
  mutate(
    correlation = map2(cell_line, phase, \(line, ph) {
      members <- strata |>
        filter(cell_line == line, phase == ph) |>
        pull(sample_id)
      # drop the phase suffix so the shared helpers see "cell_line|guide"
      subset <- strata_counts[, members]
      colnames(subset) <- str_remove(members, "\\|[^|]+$")
      subset |>
        logfc_vs_controls() |>
        correlate_targeting()
    })
  )

summarise_correlation <- function(corr) {
  values <- corr[upper.tri(corr)]
  tibble(n_guides = nrow(corr), mean_r = mean(values), min_r = min(values), max_r = max(values))
}

phase_summary <- phase_correlations |>
  mutate(stats = map(correlation, summarise_correlation)) |>
  select(-correlation) |>
  unnest(stats)

cat("\n=== guide-guide correlation within each phase ===\n")
print(phase_summary)
write_csv(phase_summary, here("output", "correlation_by_phase_summary.csv"))

# --- correlations after composition matching --------------------------------
matched_correlation <- matched_counts |>
  logfc_vs_controls() |>
  correlate_targeting()

write.csv(matched_correlation, here("output", "correlation_logfc_phase_matched.csv"))

matched_summary <- matched_correlation |>
  split_by_cell_line() |>
  imap(\(corr, line) summarise_correlation(corr) |> mutate(cell_line = line)) |>
  list_rbind() |>
  relocate(cell_line)

cat("\n=== guide-guide correlation after phase matching ===\n")
print(matched_summary)
write_csv(matched_summary, here("output", "correlation_phase_matched_summary.csv"))

# --- figures ----------------------------------------------------------------
annotation_colours <- list(
  `Cell line` = make_palette(c("A549", "H1299"), "Dark2"),
  Target = make_palette(
    c("CD274", "CD9", "FAS (CD95)", "mRuby (non-targeting)"), "Set2"
  ),
  Guide = make_palette(
    build_annotation_table(colnames(matched_counts))[["guide"]], "Paired"
  )
)

save_per_line_panels(
  matched_correlation, "heatmap_logfc_correlation_phase_matched",
  "logFC vs mRuby, phase-matched cells",
  diverging = TRUE, annotation_colours = annotation_colours
)

# publication version: titles stripped, the caption carries the description
save_per_line_panels(
  matched_correlation, "figure_guide_correlation",
  "logFC vs mRuby, phase-matched cells",
  diverging = TRUE, annotation_colours = annotation_colours, show_titles = FALSE
)

# G1-only correlations as a standalone two-panel figure, so that a figure set
# built on G1-restricted cells stays internally consistent
g1_members <- strata |> filter(phase == "G1") |> pull(sample_id)
g1_counts <- strata_counts[, g1_members]
colnames(g1_counts) <- str_remove(g1_members, "\\|[^|]+$")

g1_correlation <- g1_counts |>
  logfc_vs_controls() |>
  correlate_targeting()

write.csv(g1_correlation, here("output", "correlation_logfc_g1only.csv"))

cat("\n=== guide-guide correlation, G1 cells only ===\n")
g1_correlation |>
  split_by_cell_line() |>
  imap(\(corr, line) summarise_correlation(corr) |> mutate(cell_line = line)) |>
  list_rbind() |>
  relocate(cell_line) |>
  print()

save_per_line_panels(
  g1_correlation, "figure_guide_correlation_g1only",
  "logFC vs mRuby, G1 cells only",
  diverging = TRUE, annotation_colours = annotation_colours, show_titles = FALSE
)

# one panel per usable stratum, each on its own scale
png(here("output", "heatmap_logfc_correlation_by_phase.png"),
  width = 6.5 * nrow(phase_correlations), height = 6.5, units = "in", res = 300
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, nrow(phase_correlations))))
pwalk(phase_correlations, \(cell_line, phase, correlation, ...) {
  position <- which(
    phase_correlations[["cell_line"]] == cell_line & phase_correlations[["phase"]] == phase
  )
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = position))
  draw(
    within_line_heatmap(
      correlation, str_glue("{cell_line} - {phase} phase"),
      annotation_colours, str_glue("r ({phase})"), diverging = TRUE
    ),
    newpage = FALSE, merge_legend = TRUE,
    heatmap_legend_side = "right", annotation_legend_side = "right"
  )
  popViewport()
})
popViewport()
dev.off()
