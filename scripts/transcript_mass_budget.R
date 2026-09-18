library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)
source(here("scripts", "heatmap_helpers.R"))

# How large is the intended knockdown compared with the collateral it comes with?
# Tables only: the figure panels this once produced are not part of the final
# figure set, but the numbers are quoted in the text and so remain reproducible.
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

write_csv(budget, here("output", "transcript_mass_budget.csv"))


# --- summary quoted in the text ---------------------------------------------
summary <- budget |>
  filter(!is_control, !is.na(on_target_cpm), on_target_cpm > 0) |>
  left_join(noise_floor |> rename(rp_floor = rp_cpm), by = "cell_line") |>
  mutate(
    target_pct = 100 * on_target_cpm / 1e6,
    rp_pct = 100 * (rp_cpm - rp_floor) / 1e6,
    fold = (rp_cpm - rp_floor) / on_target_cpm
  ) |>
  select(cell_line, guide, target, target_cpm, confidence, on_target_log2,
         target_pct, rp_pct, fold) |>
  arrange(cell_line, desc(fold))

cat("\n=== transcript mass: intended target vs ribosomal collateral ===\n")
print(summary |> mutate(across(where(is.numeric), \(x) signif(x, 3))), n = Inf)
write_csv(summary, here("output", "transcript_mass_summary.csv"))
