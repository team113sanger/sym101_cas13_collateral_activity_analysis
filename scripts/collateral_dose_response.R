library(tidyverse)
library(here)
library(ggrepel)
library(patchwork)

# Does collateral activity scale with how well a guide cuts its own target?
# Uses G1-only pseudobulks so the cell cycle shift in H1299 cannot contribute.

fdr_cutoff <- 0.05

counts <- read_csv(here("output", "pseudobulk_counts_g1.csv"), show_col_types = FALSE) |>
  column_to_rownames(var = "...1") |>
  as.matrix()

logcpm <- counts |>
  prop.table(margin = 2) |>
  (\(x) log1p(x * 1e6))()

expressed <- rowMeans(exp(logcpm) - 1) >= 1
logcpm <- logcpm[expressed, ]

samples <- tibble(sample_id = colnames(logcpm)) |>
  separate_wider_delim(
    sample_id,
    delim = "|", names = c("cell_line", "guide"), cols_remove = FALSE
  ) |>
  mutate(
    family = guide |> str_to_upper() |> str_remove("_\\d+$"),
    on_target_gene = case_match(family, "CD9" ~ "CD9", "CD95" ~ "FAS", .default = NA),
    is_control = family == "MRUBY"
  )

# leave-one-out logFC against the mRuby controls of the same cell line
logfc <- samples |>
  pmap(\(sample_id, cell_line, guide, family, on_target_gene, is_control) {
    controls <- samples |>
      filter(cell_line == .env[["cell_line"]], is_control, sample_id != .env[["sample_id"]]) |>
      pull(sample_id)
    logcpm[, sample_id] - rowMeans(logcpm[, controls, drop = FALSE])
  }) |>
  set_names(samples[["sample_id"]]) |>
  bind_cols() |>
  as.matrix()
rownames(logfc) <- rownames(logcpm)

# cytoplasmic ribosomal protein genes: the module carrying the shared signal
rp_genes <- rownames(logfc) |> str_subset("^RP[LS]\\d")
cat(str_glue("ribosomal protein genes: {length(rp_genes)}"), "\n")

guide_summary <- samples |>
  mutate(
    rp_logfc = colMeans(logfc[rp_genes, sample_id]),
    collateral_magnitude = colMeans(abs(logfc[, sample_id])),
    on_target_logfc = map2_dbl(on_target_gene, sample_id, \(gene, id) {
      if (is.na(gene) || !gene %in% rownames(logfc)) NA_real_ else logfc[gene, id]
    })
  )

write_csv(guide_summary, here("output", "guide_level_collateral_summary.csv"))
print(guide_summary |> select(cell_line, guide, on_target_logfc, rp_logfc, collateral_magnitude))

# Correlation is reported with and without CD9_1, the strongest guide: it sits
# far from the others on the x axis and carries most of the leverage.
dose <- guide_summary |> filter(!is.na(on_target_logfc))

correlation_summary <- list(
  all = dose,
  without_cd9_1 = dose |> filter(str_to_upper(guide) != "CD9_1"),
  cd9_family_only = dose |> filter(family == "CD9")
) |>
  imap(\(data, subset) {
    tibble(
      subset = subset,
      n = nrow(data),
      pearson_r = cor(data[["on_target_logfc"]], data[["rp_logfc"]]),
      pearson_p = cor.test(data[["on_target_logfc"]], data[["rp_logfc"]])[["p.value"]],
      spearman_rho = cor(data[["on_target_logfc"]], data[["rp_logfc"]], method = "spearman"),
      spearman_p = cor.test(
        data[["on_target_logfc"]], data[["rp_logfc"]],
        method = "spearman", exact = FALSE
      )[["p.value"]]
    )
  }) |>
  list_rbind()

cat("\n=== on-target knockdown vs ribosomal protein depletion ===\n")
print(correlation_summary)
write_csv(correlation_summary, here("output", "dose_response_statistics.csv"))

dose_plot <- dose |>
  ggplot(aes(on_target_logfc, rp_logfc, colour = cell_line)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey40", linewidth = 0.4, alpha = 0.12) +
  geom_point(size = 2.2) +
  geom_text_repel(aes(label = guide), size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = c(A549 = "#1B9E77", H1299 = "#595959")) +
  labs(
    title = "On-target knockdown vs ribosomal protein depletion",
    subtitle = "G1 cells only; trend depends on CD9_1 (r = 0.73 with it, 0.21 without)",
    x = "on-target log2 FC (knockdown of the guide's own target)",
    y = "mean log2 FC, ribosomal protein genes",
    colour = NULL
  ) +
  theme_bw(base_size = 9)

# every gene, ranked, to show the RP module sits in the tail of the response
module_plot <- guide_summary |>
  filter(!is_control) |>
  mutate(guide_label = str_c(cell_line, guide, sep = " ")) |>
  ggplot(aes(rp_logfc, fct_reorder(guide_label, rp_logfc), fill = cell_line)) +
  geom_col() +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(A549 = "#1B9E77", H1299 = "#595959")) +
  labs(
    title = "Ribosomal protein module, per targeting guide",
    x = "mean log2 FC vs mRuby controls", y = NULL, fill = NULL
  ) +
  theme_bw(base_size = 9)

ggsave(
  here("output", "collateral_dose_response.png"),
  dose_plot + module_plot + plot_layout(widths = c(1.2, 1)),
  width = 11, height = 4.5, dpi = 300
)
