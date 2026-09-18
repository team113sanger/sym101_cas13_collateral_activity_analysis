library(tidyverse)
library(here)
library(fgsea)
library(msigdbr)
library(patchwork)
source(here("scripts", "heatmap_helpers.R"))

# Ribosomal protein depletion, guide by guide.
#
# mRuby guides are genuinely non-targeting, so they are the zero point; CD9_3
# knocks its target down by ~0.1 log2 and acts as a second, accidental inert
# control. Everything is computed on G1 cells so the phase shift in H1299
# cannot contribute.

set.seed(1)
theme_set(theme_bw(base_size = 9))
guide_type_colours <- c(
  `non-targeting (mRuby)` = "#1B9E77",
  `inert (no knockdown)` = "#7570B3",
  `targeting` = "#B2182B"
)

counts_by_phase <- read_corr(here("output", "pseudobulk_counts_by_phase.csv"))
phase_meta <- tibble(sample_id = colnames(counts_by_phase)) |>
  separate_wider_delim(
    sample_id,
    delim = "|", names = c("cell_line", "guide", "phase"), cols_remove = FALSE
  ) |>
  filter(phase == "G1") |>
  mutate(
    guide_type = case_when(
      str_detect(str_to_lower(guide), "mruby") ~ "non-targeting (mRuby)",
      str_to_upper(guide) %in% str_to_upper(ineffective_guides) ~ "inert (no knockdown)",
      .default = "targeting"
    )
  )

cpm <- prop.table(counts_by_phase[, phase_meta[["sample_id"]]], margin = 2) * 1e6
colnames(cpm) <- str_c(phase_meta[["cell_line"]], phase_meta[["guide"]], sep = "|")

# well-expressed cytoplasmic RP genes; the low-count ones only add log noise
rp_genes <- rownames(cpm) |>
  str_subset("^RP[LS]\\d") |>
  keep(\(g) mean(cpm[g, ]) >= 10)
cat(str_glue("ribosomal protein genes: {length(rp_genes)}"), "\n")

guides <- phase_meta |>
  mutate(
    sample = str_c(cell_line, guide, sep = "|"),
    rp_percent = 100 * colSums(cpm[rp_genes, sample]) / 1e6
  )

# --- panel a: absolute ribosomal content, no baseline subtracted -------------
panel_a <- guides |>
  mutate(guide = fct_reorder(str_c(cell_line, guide, sep = " "), rp_percent)) |>
  ggplot(aes(rp_percent, guide, colour = guide_type)) +
  geom_segment(aes(x = 0, xend = rp_percent, yend = guide), linewidth = 0.3, colour = "grey80") +
  geom_point(size = 2.4) +
  facet_wrap(vars(cell_line), scales = "free", ncol = 1) +
  scale_colour_manual(values = guide_type_colours, name = NULL) +
  coord_cartesian(xlim = c(20, 27.5)) +
  labs(
    title = "Absolute ribosomal protein content per guide",
    subtitle = "G1 cells, no baseline subtracted; CD9_3 shown separately as an inert control",
    x = "ribosomal protein genes, % of transcriptome", y = NULL
  )

annotation_colours <- list(
  `Cell line` = make_palette(c("A549", "H1299"), "Dark2"),
  Target = make_palette(
    c("CD274", "CD9", "FAS (CD95)", "mRuby (non-targeting)"), "Set2"
  ),
  Guide = make_palette(build_annotation_table(colnames(cpm))[["guide"]], "Paired")
)

# --- panel b: pathway-level depletion per guide ------------------------------
logfc <- guides |>
  nest(.by = cell_line) |>
  mutate(
    matrix = map2(cell_line, data, \(line, d) {
      controls <- d |>
        filter(guide_type == "non-targeting (mRuby)") |>
        pull(sample)
      d[["sample"]] |>
        set_names() |>
        map(\(s) {
          reference <- setdiff(controls, s)
          log2(cpm[, s] + 1) - rowMeans(log2(cpm[, reference, drop = FALSE] + 1))
        })
    })
  )

translation_sets <- msigdbr(
  species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"
) |>
  filter(str_detect(gs_name, "TRANSLATION|RIBOSOM|NONSENSE_MEDIATED|SRP_DEPENDENT")) |>
  (\(sets) split(sets[["gene_symbol"]], sets[["gs_name"]]))()

pathway_scores <- logfc |>
  select(cell_line, matrix) |>
  mutate(
    scores = map(matrix, \(per_guide) {
      imap(per_guide, \(ranks, sample) {
        fgsea(
          translation_sets, sort(ranks, decreasing = TRUE),
          minSize = 15, maxSize = 500, nPermSimple = 5000
        ) |>
          as_tibble() |>
          mutate(sample = sample)
      }) |>
        list_rbind()
    })
  ) |>
  select(cell_line, scores) |>
  unnest(scores) |>
  mutate(
    guide = str_split_i(sample, "\\|", 2),
    guide_type = if_else(
      str_detect(str_to_lower(guide), "mruby"), "non-targeting (mRuby)", "targeting"
    ),
    pathway = pathway |> str_remove("REACTOME_") |> str_replace_all("_", " ") |> str_trunc(40)
  )

write_csv(
  pathway_scores |> mutate(leadingEdge = map_chr(leadingEdge, \(x) str_c(x, collapse = ";"))),
  here("output", "translation_pathway_nes_per_guide.csv")
)

guide_nes <- pathway_scores |>
  summarise(
    mean_nes = mean(NES), low = min(NES), high = max(NES),
    n_depleted = sum(NES < 0 & padj < 0.05),
    .by = c(cell_line, guide, guide_type)
  )

# the RP module itself, with a confidence interval across its 104 genes, is a
# tighter read-out than averaging NES over Reactome sets that include
# mitochondrial translation and rRNA processing
module_effect <- logfc |>
  select(cell_line, matrix) |>
  mutate(
    estimate = map(matrix, \(per_guide) {
      imap(per_guide, \(ranks, sample) {
        values <- ranks[rp_genes]
        test <- t.test(values)
        tibble(
          sample = sample,
          mean_logfc = mean(values),
          conf_low = test[["conf.int"]][[1]],
          conf_high = test[["conf.int"]][[2]]
        )
      }) |>
        list_rbind()
    })
  ) |>
  select(cell_line, estimate) |>
  unnest(estimate) |>
  mutate(
    guide = str_split_i(sample, "\\|", 2),
    guide_type = case_when(
      str_detect(str_to_lower(guide), "mruby") ~ "non-targeting (mRuby)",
      str_to_upper(guide) %in% str_to_upper(ineffective_guides) ~ "inert (no knockdown)",
      .default = "targeting"
    )
  )

write_csv(module_effect, here("output", "rp_module_effect_per_guide.csv"))

panel_b <- module_effect |>
  ggplot(aes(mean_logfc, fct_reorder(guide, mean_logfc), colour = guide_type)) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.3) +
  geom_errorbarh(aes(xmin = conf_low, xmax = conf_high), height = 0, linewidth = 0.4) +
  geom_point(size = 2.4) +
  facet_wrap(vars(cell_line), scales = "free_y", ncol = 1) +
  scale_colour_manual(values = guide_type_colours, name = NULL) +
  labs(
    title = "Ribosomal protein module shift per guide",
    subtitle = str_glue("mean log2 FC over {length(rp_genes)} genes, 95% CI"),
    x = "mean log2 FC vs mRuby controls", y = NULL
  )

write_csv(guide_nes, here("output", "translation_nes_summary_per_guide.csv"))

# pathway x guide detail
nes_heatmaps <- pathway_scores |>
  nest(.by = cell_line) |>
  mutate(
    heatmap = map2(cell_line, data, \(line, d) {
      values <- d |>
        select(pathway, guide, NES) |>
        pivot_wider(names_from = guide, values_from = NES) |>
        column_to_rownames("pathway") |>
        as.matrix()
      limit <- max(abs(values), na.rm = TRUE)
      Heatmap(
        values,
        name = str_glue("NES ({line})"),
        col = colorRamp2(c(-limit, 0, limit), c("#2166AC", "white", "#B2182B")),
        column_title = str_glue("{line}"),
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 7),
        top_annotation = build_annotation(
          build_annotation_table(str_c(line, unique(d[["guide"]]), sep = "|")),
          "column", annotation_colours, fields = c("Target", "Guide")
        ),
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
      )
    })
  )

ggsave(
  here("output", "rp_depletion_by_guide.png"),
  panel_a + panel_b + plot_layout(guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(legend.position = "bottom"),
  width = 12, height = 7.5, dpi = 300
)

# --- gene-level view: every RP gene, every guide -----------------------------
rp_heatmaps <- logfc |>
  mutate(
    heatmap = map2(cell_line, matrix, \(line, per_guide) {
      values <- per_guide |>
        map(\(x) x[rp_genes]) |>
        bind_cols() |>
        as.matrix()
      rownames(values) <- rp_genes
      # CD9_1 is far stronger than the rest; cap at the 98th percentile so the
      # other guides are not washed out (values beyond the cap saturate)
      limit <- quantile(abs(values), 0.98)

      Heatmap(
        values,
        name = str_glue("log2 FC ({line})"),
        col = colorRamp2(c(-limit, 0, limit), c("#2166AC", "white", "#B2182B")),
        column_title = str_glue("{line} - {length(rp_genes)} ribosomal protein genes"),
        column_title_gp = gpar(fontsize = 11, fontface = "bold"),
        top_annotation = build_annotation(
          build_annotation_table(colnames(values)), "column", annotation_colours,
          fields = c("Target", "Guide")
        ),
        show_row_names = FALSE,
        column_names_gp = gpar(fontsize = 7),
        column_labels = str_split_i(colnames(values), "\\|", 2),
        show_row_dend = FALSE,
        heatmap_legend_param = list(title_gp = gpar(fontsize = 8), labels_gp = gpar(fontsize = 7))
      )
    })
  )

png(here("output", "rp_gene_heatmap_by_guide.png"),
  width = 13, height = 7, units = "in", res = 300
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, nrow(rp_heatmaps))))
walk(seq_len(nrow(rp_heatmaps)), \(i) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
  draw(
    rp_heatmaps[["heatmap"]][[i]],
    newpage = FALSE, merge_legend = TRUE,
    heatmap_legend_side = "right", annotation_legend_side = "right"
  )
  popViewport()
})
popViewport()
dev.off()

png(here("output", "rp_pathway_nes_heatmap.png"),
  width = 15, height = 6, units = "in", res = 300
)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, nrow(nes_heatmaps))))
walk(seq_len(nrow(nes_heatmaps)), \(i) {
  pushViewport(viewport(layout.pos.row = 1, layout.pos.col = i))
  draw(
    nes_heatmaps[["heatmap"]][[i]],
    newpage = FALSE, merge_legend = TRUE,
    heatmap_legend_side = "right", annotation_legend_side = "right"
  )
  popViewport()
})
popViewport()
dev.off()
