library(tidyverse)
library(here)
library(ggrepel)

# Pathway responses computed independently in each cell line, plotted against
# each other. Nothing forces a point onto the diagonal, so proximity to it is a
# replication result rather than a property of the plot.

theme_set(theme_bw(base_size = 9))

gsea <- read_csv(here("output", "gsea_targeting_vs_control.csv"), show_col_types = FALSE)

pathways <- gsea |>
  filter(contrast %in% c("A549_g1only", "H1299_g1only"), collection == "reactome") |>
  select(contrast, pathway, NES, padj) |>
  pivot_wider(names_from = contrast, values_from = c(NES, padj)) |>
  drop_na() |>
  mutate(
    both_significant = padj_A549_g1only < 0.05 & padj_H1299_g1only < 0.05,
    translation = str_detect(pathway, "TRANSLATION|RIBOSOM|NONSENSE_MEDIATED|SRP_DEPENDENT")
  )

significant <- pathways |> filter(both_significant)
cat(str_glue(
  "sets plotted: {nrow(pathways)} | significant in both: {nrow(significant)} | ",
  "concordant: {sum(sign(significant$NES_A549_g1only) == sign(significant$NES_H1299_g1only))} | ",
  "depleted in both: {sum(significant$NES_A549_g1only < 0 & significant$NES_H1299_g1only < 0)}"
), "\n")
cat(str_glue(
  "NES r among significant: {round(cor(significant$NES_A549_g1only, significant$NES_H1299_g1only), 3)} | ",
  "across all sets: {round(cor(pathways$NES_A549_g1only, pathways$NES_H1299_g1only), 3)}"
), "\n")

scatter <- pathways |>
  ggplot(aes(NES_A549_g1only, NES_H1299_g1only)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_vline(xintercept = 0, linewidth = 0.3) +
  geom_abline(linetype = "dashed", linewidth = 0.3, colour = "grey50") +
  geom_point(data = \(d) filter(d, !both_significant), colour = "grey85", size = 0.5) +
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
  labs(x = "NES (A549)", y = "NES (H1299)")

ggsave(here("output", "figure_nes_concordance.png"), scatter,
  width = 7.2, height = 5.6, dpi = 300)
ggsave(here("output", "figure_nes_concordance.pdf"), scatter,
  width = 7.2, height = 5.6, device = cairo_pdf)
