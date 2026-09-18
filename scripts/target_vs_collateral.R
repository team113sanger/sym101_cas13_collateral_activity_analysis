library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
source(here("scripts", "heatmap_helpers.R"))

# How large is the intended knockdown compared with the collateral it comes with?
#
# Everything is expressed as parts per million of the transcriptome, so the two
# can be added up on one scale. mRuby controls are carried through the identical
# calculation to give the noise floor, since a pseudobulk of a few hundred cells
# shows thousands of "shifted" genes by chance alone.
#
# Guides are tiered by how well their target can be measured in that cell line:
#   high  - target >= 10 CPM (CD9 in both lines, FAS in A549)
#   low   - target 1-10 CPM (FAS in H1299, 3.7 CPM)
#   none  - target < 1 CPM (CD274, ~0.8 CPM, 1-9 counts per pseudobulk)
# CD274 is below the expression filter used for the pseudobulk matrices, so its
# counts are read back from the h5ad separately.

theme_set(theme_bw(base_size = 9))

counts_by_phase <- read_corr(here("output", "pseudobulk_counts_by_phase.csv"))
phase_meta <- tibble(sample_id = colnames(counts_by_phase)) |>
  separate_wider_delim(
    sample_id,
    delim = "|", names = c("cell_line", "guide", "phase"), cols_remove = FALSE
  ) |>
  filter(phase == "G1")

cpm <- prop.table(counts_by_phase[, phase_meta[["sample_id"]]], margin = 2) * 1e6
colnames(cpm) <- str_c(phase_meta[["cell_line"]], phase_meta[["guide"]], sep = "|")

rp_genes <- rownames(cpm) |>
  str_subset("^RP[LS]\\d") |>
  keep(\(g) mean(cpm[g, ]) >= 10)

# CD274 sits below the expression filter, so its per-pseudobulk CPM comes from
# the raw counts; the low totals are kept so the figure can flag them
cd274 <- read_csv(here("output", "cd274_counts_g1.csv"), show_col_types = FALSE) |>
  mutate(sample = str_c(cell_line, guide, sep = "|"))

target_cpm_in_line <- function(gene, line) {
  columns <- colnames(cpm) |> str_subset(str_glue("^{line}\\|"))
  if (gene %in% rownames(cpm)) mean(cpm[gene, columns]) else 0
}

shift_threshold <- 0.1

budget <- phase_meta |>
  distinct(cell_line) |>
  mutate(
    rows = map(cell_line, \(line) {
      samples <- colnames(cpm) |> str_subset(str_glue("^{line}\\|"))
      controls <- samples |> str_subset(regex("mruby", ignore_case = TRUE))

      samples |>
        discard(is_ineffective) |>
        map(\(sample) {
          reference <- setdiff(controls, sample)
          baseline <- rowMeans(cpm[, reference, drop = FALSE])
          observed <- cpm[, sample]
          delta <- observed - baseline
          log_change <- log2(observed + 1) - log2(baseline + 1)

          guide <- str_split_i(sample, "\\|", 2)
          target <- guide |>
            str_to_upper() |>
            str_remove("_\\d+$") |>
            recode(CD95 = "FAS")

          # measured within this cell line, never pooled across lines
          target_level <- target_cpm_in_line(target, line)
          if (target == "CD274") {
            target_level <- mean(cd274[["cd274_cpm"]][cd274[["cell_line"]] == line])
          }

          on_target_log2 <- NA_real_
          on_target_cpm <- NA_real_
          if (target %in% rownames(cpm)) {
            on_target_log2 <- log_change[[target]]
            on_target_cpm <- -delta[[target]]
          } else if (target == "CD274") {
            observed_cd274 <- cd274[["cd274_cpm"]][cd274[["sample"]] == sample]
            reference_cd274 <- cd274[["cd274_cpm"]][cd274[["sample"]] %in% reference]
            on_target_log2 <- log2(observed_cd274 + 1) - log2(mean(reference_cd274) + 1)
            on_target_cpm <- mean(reference_cd274) - observed_cd274
          }

          tibble(
            guide = guide,
            is_control = str_detect(str_to_lower(guide), "mruby"),
            target = target,
            target_cpm = target_level,
            confidence = case_when(
              target_level >= 10 ~ "target well measured",
              target_level >= 1 ~ "target poorly measured",
              .default = "target not quantifiable"
            ),
            on_target_log2 = on_target_log2,
            on_target_cpm = on_target_cpm,
            rp_cpm = -sum(delta[rp_genes]),
            redistributed_cpm = sum(abs(delta)) / 2,
            genes_shifted = sum(abs(log_change) > shift_threshold)
          )
        }) |>
        list_rbind()
    })
  ) |>
  unnest(rows)

noise_floor <- budget |>
  filter(is_control) |>
  summarise(
    rp_cpm = mean(rp_cpm),
    redistributed_cpm = mean(redistributed_cpm),
    genes_shifted = mean(genes_shifted),
    .by = cell_line
  )

cat("=== control noise floor (mRuby vs the other mRuby guides) ===\n")
print(noise_floor)

write_csv(budget, here("output", "target_vs_collateral_budget.csv"))

# --- panel a: transcript mass, target vs collateral --------------------------
mass <- budget |>
  filter(!is_control, !is.na(on_target_cpm)) |>
  left_join(noise_floor |> rename(rp_floor = rp_cpm, redist_floor = redistributed_cpm), by = "cell_line") |>
  mutate(
    `intended target` = on_target_cpm,
    `ribosomal module` = rp_cpm - rp_floor,
    label = str_c(cell_line, guide, sep = " ")
  ) |>
  select(cell_line, guide, label, confidence, on_target_log2, `intended target`, `ribosomal module`) |>
  pivot_longer(
    c(`intended target`, `ribosomal module`),
    names_to = "component", values_to = "cpm"
  )

panel_a <- mass |>
  filter(cpm > 0) |>
  mutate(label = str_c(label, if_else(confidence == "target well measured", "", " *"))) |>
  ggplot(aes(cpm, fct_reorder(label, cpm), fill = component)) +
  geom_col(
    aes(alpha = confidence),
    position = position_dodge(width = 0.7), width = 0.65
  ) +
  scale_alpha_manual(
    values = c(
      `target well measured` = 1,
      `target poorly measured` = 0.55,
      `target not quantifiable` = 0.25
    ), name = NULL
  ) +
  scale_x_log10(
    labels = scales::label_comma(),
    sec.axis = sec_axis(~ . / 1e4, name = "% of transcriptome")
  ) +
  scale_fill_manual(values = c(`intended target` = "#4292C6", `ribosomal module` = "#B2182B"), name = NULL) +
  labs(
    title = "Transcript mass removed: intended target vs ribosomal collateral",
    subtitle = "G1 cells; faded bars (*) have unquantifiable targets",
    x = "parts per million of the transcriptome (log scale)", y = NULL
  )

# --- panel b: the ratio ------------------------------------------------------
ratio <- budget |>
  filter(!is_control, !is.na(on_target_cpm), on_target_cpm > 0) |>
  left_join(noise_floor |> rename(rp_floor = rp_cpm), by = "cell_line") |>
  mutate(fold = (rp_cpm - rp_floor) / on_target_cpm, label = str_c(cell_line, guide, sep = " "))

panel_b <- ratio |>
  filter(fold > 0) |>
  mutate(label = str_c(label, " (", target, " ", round(target_cpm, 1), " CPM)")) |>
  ggplot(aes(fold, fct_reorder(label, fold), colour = cell_line)) +
  geom_segment(aes(x = 1, xend = fold, yend = fct_reorder(label, fold)), linewidth = 0.3, colour = "grey75") +
  geom_point(aes(shape = confidence), size = 2.6) +
  scale_shape_manual(
    values = c(
      `target well measured` = 16,
      `target poorly measured` = 1,
      `target not quantifiable` = 4
    ), name = NULL
  ) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.3) +
  scale_x_log10() +
  scale_colour_manual(values = c(A549 = "#1B9E77", H1299 = "#7570B3"), name = NULL) +
  labs(
    title = "Collateral costs more transcript than the knockdown achieves",
    subtitle = "RP loss / on-target loss; dashed line = parity",
    x = "fold more transcript lost from the ribosomal module (log scale)", y = NULL
  )

# --- panel c: where the target sits in the whole distribution ----------------
strongest <- c("A549|CD9_1", "H1299|CD9_1")

# Genes are re-filtered on the G1 columns actually plotted, matching
# logfc_vs_controls() in heatmap_helpers.R. Without this the panel inherits the
# filter computed on all cells and shows ~195 S/G2M-specific transcripts that
# fall below 1 CPM in G1 and carry unstable log ratios.
expressed_g1 <- rownames(cpm)[rowMeans(cpm) >= 1]
cat("\ngenes passing mean CPM >= 1 in G1 cells:", length(expressed_g1), "\n")

ranked <- strongest |>
  set_names() |>
  map(\(sample) {
    line <- str_split_i(sample, "\\|", 1)
    controls <- colnames(cpm) |>
      str_subset(str_glue("^{line}\\|")) |>
      str_subset(regex("mruby", ignore_case = TRUE))
    baseline <- rowMeans(cpm[expressed_g1, controls, drop = FALSE])
    log_change <- log2(cpm[expressed_g1, sample] + 1) - log2(baseline + 1)

    tibble(
      sample = sample,
      cell_line = line,
      gene = names(log_change),
      log2fc = log_change,
      class = case_when(
        gene == "CD9" ~ "CD9 (intended target)",
        gene %in% rp_genes ~ "ribosomal protein",
        .default = "all other genes"
      )
    )
  }) |>
  list_rbind()

panel_c <- ranked |>
  ggplot(aes(log2fc, fill = class)) +
  geom_histogram(
    data = \(d) filter(d, class == "all other genes"),
    bins = 120, alpha = 0.55
  ) +
  geom_rug(
    data = \(d) filter(d, class == "ribosomal protein"),
    aes(colour = class), sides = "b", linewidth = 0.3, length = unit(0.05, "npc")
  ) +
  geom_vline(
    data = \(d) filter(d, class == "CD9 (intended target)"),
    aes(xintercept = log2fc), colour = "#4292C6", linewidth = 0.7
  ) +
  geom_text(
    data = \(d) filter(d, class == "CD9 (intended target)"),
    aes(x = log2fc, y = Inf, label = str_glue("CD9  {round(log2fc, 1)}")),
    colour = "#4292C6", hjust = -0.1, vjust = 1.8, size = 2.6
  ) +
  facet_wrap(vars(cell_line), scales = "free_y") +
  coord_cartesian(xlim = c(-3.2, 1.5)) +
  scale_fill_manual(values = c(`all other genes` = "grey70"), guide = "none") +
  scale_colour_manual(values = c(`ribosomal protein` = "#B2182B"), name = NULL) +
  labs(
    title = "The target is the single largest change; the collateral is the wider one",
    subtitle = "CD9_1, per-gene log2 FC vs mRuby controls. Red ticks mark ribosomal protein genes",
    x = "log2 FC", y = "genes"
  )

ggsave(
  here("output", "target_vs_collateral.png"),
  (panel_a + panel_b) / panel_c + plot_annotation(tag_levels = "a"),
  width = 13, height = 8.5, dpi = 300
)

# standalone version of the distribution panel, titles stripped so the caption
# carries the description instead
histogram_panel <- panel_c + labs(title = NULL, subtitle = NULL)
ggsave(
  here("output", "figure_target_vs_rp_histogram.png"),
  histogram_panel, width = 8, height = 3.6, dpi = 300
)
# vector version for figure assembly; cairo_pdf embeds the fonts
ggsave(
  here("output", "figure_target_vs_rp_histogram.pdf"),
  histogram_panel, width = 8, height = 3.6, device = cairo_pdf
)

cat("\n=== panel c contents ===\n")
cat("genes plotted:", n_distinct(ranked[["gene"]]), "\n")
cat("ribosomal protein genes marked:", length(rp_genes), "\n")
print(ranked |> filter(class == "CD9 (intended target)") |> select(cell_line, gene, log2fc))

cat("\n=== headline numbers (CD9 guides, G1) ===\n")
print(
  ratio |>
    mutate(
      target_pct = 100 * on_target_cpm / 1e6,
      rp_pct = 100 * (rp_cpm - rp_floor) / 1e6
    ) |>
    select(cell_line, guide, target, target_cpm, confidence, on_target_log2, target_pct, rp_pct, fold) |>
    mutate(across(where(is.numeric), \(x) round(x, 3)))
)
