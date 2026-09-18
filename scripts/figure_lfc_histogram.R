library(tidyverse)
library(here)
source(here("scripts", "heatmap_helpers.R"))

# Where the intended knockdown sits in the distribution of all per-gene changes,
# and where the ribosomal protein genes sit relative to it. CD9_1 is used because
# CD9 is the only target expressed highly enough to quantify reliably.

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

# re-filtered on the G1 columns actually plotted, matching logfc_vs_controls();
# without this the panel inherits the filter computed on all cells and shows
# ~195 S/G2M-specific transcripts that fall below 1 CPM in G1
expressed <- rownames(cpm)[rowMeans(cpm) >= 1]
rp_genes <- expressed |>
  str_subset("^RP[LS]\\d") |>
  keep(\(g) mean(cpm[g, ]) >= 10)

cat(str_glue("genes plotted: {length(expressed)} | ribosomal protein genes: {length(rp_genes)}"), "\n")

ranked <- c("A549|CD9_1", "H1299|CD9_1") |>
  set_names() |>
  map(\(sample) {
    line <- str_split_i(sample, "\\|", 1)
    controls <- colnames(cpm) |>
      str_subset(str_glue("^{line}\\|")) |>
      str_subset(regex("mruby", ignore_case = TRUE))
    baseline <- rowMeans(cpm[expressed, controls, drop = FALSE])
    log_change <- log2(cpm[expressed, sample] + 1) - log2(baseline + 1)

    tibble(
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

print(ranked |> filter(class == "CD9 (intended target)") |> select(cell_line, gene, log2fc))

histogram <- ranked |>
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
  labs(x = "log2 FC", y = "genes")

ggsave(here("output", "figure_target_vs_rp_histogram.png"), histogram,
  width = 8, height = 3.6, dpi = 300)
ggsave(here("output", "figure_target_vs_rp_histogram.pdf"), histogram,
  width = 8, height = 3.6, device = cairo_pdf)
