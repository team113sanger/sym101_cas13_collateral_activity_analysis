library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
source(here("scripts", "heatmap_helpers.R"))

# Figures for the three-part collateral story:
#   1. a target-independent transcriptional change
#   2. that change is a translation/ribosome program
#   3. how its severity compares between the two lines
#
# A549 and H1299 are always compared in G1, the only phase both lines have in
# useful numbers, so proliferation state is not driving the contrast.

theme_set(theme_bw(base_size = 9))
line_colours <- c(A549 = "#1B9E77", H1299 = "#7570B3")

counts_by_phase <- read_corr(here("output", "pseudobulk_counts_by_phase.csv"))
phase_meta <- tibble(sample_id = colnames(counts_by_phase)) |>
  separate_wider_delim(
    sample_id,
    delim = "|", names = c("cell_line", "guide", "phase"), cols_remove = FALSE
  )

logfc_for_stratum <- function(line, target_phase) {
  members <- phase_meta |>
    filter(cell_line == line, phase == target_phase) |>
    pull(sample_id)
  subset <- counts_by_phase[, members]
  colnames(subset) <- str_remove(members, "\\|[^|]+$")
  logfc_vs_controls(subset)
}

g1_logfc <- c("A549", "H1299") |>
  set_names() |>
  map(\(line) logfc_for_stratum(line, "G1"))

# ---------------------------------------------------------------------------
# Figure 1 - the change is target-independent
# ---------------------------------------------------------------------------

guide_target <- function(ids) {
  ids |>
    str_split_i("\\|", 2) |>
    str_to_upper() |>
    str_remove("_\\d+$") |>
    recode(CD95 = "FAS")
}

pair_table <- g1_logfc |>
  imap(\(logfc, line) {
    targeting <- colnames(logfc) |>
      str_subset("mruby", negate = TRUE) |>
      discard(is_ineffective)
    corr <- cor(logfc[, targeting])
    which(upper.tri(corr), arr.ind = TRUE) |>
      as_tibble() |>
      mutate(
        cell_line = line,
        guide_a = colnames(corr)[row],
        guide_b = colnames(corr)[col],
        r = corr[cbind(row, col)],
        pair_type = if_else(
          guide_target(guide_a) == guide_target(guide_b),
          "same target", "different target"
        )
      )
  }) |>
  list_rbind()

pair_tests <- pair_table |>
  nest(.by = cell_line) |>
  mutate(
    p = map_dbl(data, \(d) wilcox.test(r ~ pair_type, data = d)[["p.value"]]),
    label = str_glue("Wilcoxon p = {round(p, 2)}")
  )

panel_1a <- pair_table |>
  ggplot(aes(pair_type, r, colour = cell_line)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "grey40") +
  geom_jitter(width = 0.12, size = 1.4, alpha = 0.85) +
  geom_text(
    data = pair_tests, aes(x = 1.5, y = Inf, label = label),
    inherit.aes = FALSE, vjust = 1.6, size = 2.8
  ) +
  facet_wrap(vars(cell_line)) +
  scale_colour_manual(values = line_colours, guide = "none") +
  labs(
    title = "Guides against the same gene are no more alike than guides against different genes",
    subtitle = "pairwise correlation of guide logFC profiles, G1 cells",
    x = NULL, y = "Pearson r between guide pairs"
  )

# on-target knockdown is unrelated to how much collateral a guide produces
knockdown <- g1_logfc |>
  imap(\(logfc, line) {
    targeting <- colnames(logfc) |> str_subset("mruby", negate = TRUE)
    rp_genes <- rownames(logfc) |> str_subset("^RP[LS]\\d")
    tibble(
      cell_line = line,
      guide = str_split_i(targeting, "\\|", 2),
      target = guide_target(targeting),
      on_target = map2_dbl(targeting, target, \(id, gene) {
        if (gene %in% rownames(logfc)) logfc[gene, id] else NA_real_
      }),
      rp_logfc = colMeans(logfc[rp_genes, targeting])
    )
  }) |>
  list_rbind()

panel_1b <- knockdown |>
  drop_na(on_target) |>
  ggplot(aes(on_target, rp_logfc, colour = cell_line)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_point(size = 2.2) +
  geom_text_repel(aes(label = guide), size = 2.5, show.legend = FALSE) +
  scale_colour_manual(values = line_colours) +
  labs(
    title = "Severity tracks guide potency, not which gene is targeted",
    subtitle = str_c(
      "CD9_1 is the strongest guide and most depleted in both lines; CD9_3 cuts ",
      "nothing and shows none. Trend rests heavily on CD9_1 (r = 0.21 without it)."
    ),
    x = "on-target log2 FC", y = "ribosomal protein module log2 FC", colour = NULL
  )

ggsave(
  here("output", "figure1_target_independence.png"),
  panel_1a / panel_1b + plot_annotation(tag_levels = "a"),
  width = 9, height = 8, dpi = 300
)

# ---------------------------------------------------------------------------
# Figure 2 - the shared change is a translation / ribosome program
# ---------------------------------------------------------------------------

gsea <- read_csv(here("output", "gsea_targeting_vs_control.csv"), show_col_types = FALSE)

shared_pathways <- gsea |>
  filter(contrast %in% c("A549_g1only", "H1299_g1only"), collection == "reactome") |>
  select(contrast, pathway, NES, padj) |>
  pivot_wider(names_from = contrast, values_from = c(NES, padj)) |>
  drop_na() |>
  mutate(
    both_significant = padj_A549_g1only < 0.05 & padj_H1299_g1only < 0.05,
    translation = str_detect(pathway, "TRANSLATION|RIBOSOM|NONSENSE_MEDIATED|SRP_DEPENDENT")
  )

panel_2a <- shared_pathways |>
  ggplot(aes(NES_A549_g1only, NES_H1299_g1only)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(
    data = \(d) filter(d, !both_significant),
    colour = "grey85", size = 0.5
  ) +
  geom_point(
    data = \(d) filter(d, both_significant),
    aes(colour = translation), size = 2.1
  ) +
  geom_text_repel(
    data = \(d) filter(d, translation, both_significant) |> slice_min(NES_A549_g1only, n = 9),
    aes(label = str_remove(pathway, "REACTOME_") |> str_replace_all("_", " ") |> str_trunc(34)),
    size = 2.1, seed = 1, max.overlaps = Inf, min.segment.length = 0,
    segment.size = 0.2, segment.colour = "grey55",
    box.padding = 0.42, point.padding = 0.25,
    nudge_x = 2.1, hjust = 0, direction = "y", xlim = c(-0.4, NA)
  ) +
  scale_colour_manual(
    values = c(`TRUE` = "#B2182B", `FALSE` = "grey75"),
    labels = c(`TRUE` = "translation / ribosome", `FALSE` = "other"), name = NULL
  ) +
  labs(
    title = "Where both lines respond, they respond the same way",
    subtitle = str_c(
      "24 of 1015 Reactome sets are significant in both lines; all 24 agree in ",
      "direction (NES r = 0.83). Translation sets are the strongest of them. ",
      "Grey = not significant in both."
    ),
    x = "NES (A549)", y = "NES (H1299)"
  )

# where the ribosomal protein genes sit in the overall response
gene_distribution <- g1_logfc |>
  imap(\(logfc, line) {
    targeting <- colnames(logfc) |> str_subset("mruby", negate = TRUE)
    tibble(
      cell_line = line,
      gene = rownames(logfc),
      mean_logfc = rowMeans(logfc[, discard(targeting, is_ineffective)])
    ) |>
      mutate(module = if_else(str_detect(gene, "^RP[LS]\\d"), "ribosomal protein", "all other genes"))
  }) |>
  list_rbind()

panel_2b <- gene_distribution |>
  ggplot(aes(mean_logfc, colour = module)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  stat_ecdf(linewidth = 0.7) +
  coord_cartesian(xlim = c(-0.25, 0.25)) +
  facet_wrap(vars(cell_line)) +
  scale_colour_manual(values = c(`ribosomal protein` = "#B2182B", `all other genes` = "grey55"), name = NULL) +
  labs(
    title = "Ribosomal protein genes shift down relative to the transcriptome",
    subtitle = "cumulative distribution of mean log2 FC across targeting guides, G1 cells",
    x = "mean log2 FC (targeting vs mRuby)", y = "cumulative fraction of genes"
  )

ggsave(
  here("output", "figure2_translation_program.png"),
  panel_2a / panel_2b + plot_annotation(tag_levels = "a"),
  width = 9, height = 8.5, dpi = 300
)

# standalone NES scatter, titles stripped so the caption carries the description
nes_panel <- panel_2a + labs(title = NULL, subtitle = NULL)
ggsave(
  here("output", "figure_nes_concordance.png"),
  nes_panel, width = 7.2, height = 5.6, dpi = 300
)
ggsave(
  here("output", "figure_nes_concordance.pdf"),
  nes_panel, width = 7.2, height = 5.6, device = cairo_pdf
)

# ---------------------------------------------------------------------------
# Figure 3 - severity between lines, with the proliferation confound shown
# ---------------------------------------------------------------------------

severity <- read_csv(here("output", "rp_module_severity.csv"), show_col_types = FALSE)

panel_3a <- severity |>
  filter(phase == "G1") |>
  ggplot(aes(group, rp_total_pct, colour = cell_line)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_boxplot(outlier.shape = NA, width = 0.5, colour = "grey40") +
  geom_jitter(width = 0.12, size = 1.8) +
  facet_wrap(vars(cell_line)) +
  scale_colour_manual(values = line_colours, guide = "none") +
  labs(
    title = "Severity in G1: directionally milder in A549, but not significantly",
    subtitle = str_c(
      "median -1.5% (A549) vs -2.8% (H1299), Wilcoxon p = 0.51 on 8 guides per ",
      "line. Distributions overlap; CD9_1 reaches -9% in both."
    ),
    x = NULL, y = "ribosomal protein module change (%)"
  )

phase_summary <- read_csv(here("output", "correlation_by_phase_summary.csv"), show_col_types = FALSE)

panel_3b <- phase_summary |>
  mutate(stratum = str_c(cell_line, phase, sep = " ")) |>
  ggplot(aes(fct_reorder(stratum, mean_r), mean_r, fill = cell_line)) +
  geom_col(width = 0.6) +
  geom_errorbar(aes(ymin = min_r, ymax = max_r), width = 0.15, linewidth = 0.3) +
  coord_flip() +
  scale_fill_manual(values = line_colours, name = NULL) +
  labs(
    title = "Consistency rises with proliferation, not cell line",
    subtitle = str_c(
      "in G1, A549 (0.27) and H1299 (0.26) are equal; H1299 only pulls ahead in ",
      "S and G2M, which A549 barely has"
    ),
    x = NULL, y = "mean guide-guide correlation (bars show range)"
  )

# the confound itself
phase_mix <- read_csv(here("output", "single_cell_dose.csv"), show_col_types = FALSE) |>
  count(cell_line, phase) |>
  mutate(fraction = n / sum(n), .by = cell_line) |>
  mutate(phase = fct_relevel(phase, "G1", "S", "G2M"))

panel_3c <- phase_mix |>
  ggplot(aes(cell_line, fraction, fill = phase)) +
  geom_col(width = 0.6) +
  scale_fill_brewer(palette = "Blues", name = NULL) +
  scale_y_continuous(labels = scales::percent) +
  labs(
    title = "Why pooled severity comparisons mislead",
    subtitle = "A549 is 91% G1; H1299 is mostly cycling",
    x = NULL, y = NULL
  )

ggsave(
  here("output", "figure3_severity_between_lines.png"),
  panel_3a / (panel_3b + panel_3c + plot_layout(widths = c(1.6, 1))) +
    plot_annotation(tag_levels = "a"),
  width = 10, height = 8, dpi = 300
)
